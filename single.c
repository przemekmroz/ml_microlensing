#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "asa047.h"
#include "find_transient.h"
#include "single.h"

#define BLENDING_LIMIT -0.1
#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

/* Identifying microlensing events using deep neural nets
 * 
 * Routines for fitting PSPL models.
 * 
 * P. Mroz @ Caltech, 1 Aug 2020
 */

double magnification (double pars[], double t) {
    
    /* PSPL magnification: t0 = pars[0]; tE = pars[1]; u0 = pars[2] */
    double tau,u,A;
    
    tau = (t-pars[0])/pars[1];
    u = sqrt(tau*tau + pars[2]*pars[2]);
    A = (u*u+2.0)/sqrt(u*u+4.0)/u;
    
    return A;
}

void get_blending (double *params, double *time, double *flux, double 
    *f_err, int *outlier, int ndata, double *Fs, double *Fb, double 
    *resid, double *chi2) {
    
    /* Calculating blending parameters given data */
    
    double A[ndata],w,W0,W1,W2;
    double sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0,sum5=0.0;
    int i;
        
    for (i=0; i<ndata; i++) {
        A[i] = magnification(params,time[i]);
        if (outlier[i] == 1) continue;
        w = 1.0/(f_err[i]*f_err[i]);
        sum1 += A[i]*A[i]*w;
        sum2 += A[i]*w;
        sum3 += w;
        sum4 += flux[i]*A[i]*w;
        sum5 += flux[i]*w;
    }
    W0 = sum1*sum3-sum2*sum2;
    W1 = sum4*sum3-sum2*sum5;
    W2 = sum1*sum5-sum4*sum2;
    
    if (W0 != 0.0) {
        *Fs = W1/W0;
        *Fb = W2/W0;
    } else {
        *Fs = -1.0;
        *Fb = 0.0;
    }
    
    *chi2 = 0.0;
    for (i=0; i<ndata; i++) {
        w = flux[i]-(*Fs)*A[i]-(*Fb);
        w /= f_err[i];
        w = w*w;
        resid[i] = w;
        if (outlier[i] == 0) *chi2 += w;
    }
    
}

void get_blending_fixed (double *params, double *time, double *flux, 
    double *f_err, int *outlier, int ndata, double *Fs, double *Fb, 
    double *resid, double *chi2) {
    
    /* Calculating blending parameters given data (assuming Fb=0) */
    
    double A[ndata],w;
    double sum1=0.0,sum2=0.0;
    int i;
        
    for (i=0; i<ndata; i++) {
        A[i] = magnification(params,time[i]);
        if (outlier[i] == 1) continue;
        w = 1.0/(f_err[i]*f_err[i]);
        sum1 += A[i]*A[i]*w;
        sum2 += flux[i]*A[i]*w;
    }
    
    *Fs = sum2/sum1;
    *Fb = 0.0;
    
    *chi2 = 0.0;
    for (i=0; i<ndata; i++) {
        w = flux[i]-(*Fs)*A[i]-(*Fb);
        w /= f_err[i];
        w = w*w;
        resid[i] = w;
        if (outlier[i] == 0) *chi2 += w;
    }
    
}

double chi2 (double *params, double *time, double *flux, double *f_err, 
    int *outlier, int ndata) {
    
    /* Calculating chi^2 */
    
    double resid[ndata],Fs,Fb,res;
    
    get_blending(params,time,flux,f_err,outlier,ndata,&Fs,&Fb,resid,&res);
    
    return res;
    
}

double chi2_fixed (double *params, double *time, double *flux, double 
    *f_err, int *outlier, int ndata) {
    
    /* Calculating chi^2 */
    
    double resid[ndata],Fs,Fb,res;
    
    get_blending_fixed(params,time,flux,f_err,outlier,ndata,&Fs,&Fb,resid,&res);
    
    return res;
    
}

int remove_4_sigma_outliers (double *time, double *resid, int *flag, 
    int ndata) {
    
    /* Flagging outliers in residuals from the PSPL fit */
    
    int i,nbad=0;
    
    if (resid[0] > 16.0 && resid[1] < 1.0 && flag[0] == 0) {
        flag[0] = 1; 
        nbad += 1;
    }
    for (i=1; i<(ndata-1); i++) {
        if (resid[i] > 16.0 && 
            (resid[i-1]+resid[i+1]) < 9.0 && 
            flag[i] == 0 && 
            ((time[i]-time[i-1]) <= 1.5 || (time[i+1]-time[i]) <= 1.5) 
           ) {
               flag[i] = 1; 
               nbad += 1;
             }
    }
    if (resid[ndata-1] > 16.0 && resid[ndata-2] < 1.0 && flag[ndata-1] == 0) {
        flag[ndata-1] = 1; 
        nbad += 1;
    }
 
    return nbad;
}

void calculate_chi2 (double *chi2, int *ndata, double *time, double 
    *resid, int *outlier, int npts, double tmin, double tmax) {
    
    /* Calculate chi^2 for data points in the time range [tmin,tmax] */
    
    int i;

    *chi2 = 0.0;
    *ndata = 0;

    for (i=0; i<npts; i++) {
        if (outlier[i] == 1) continue;
        if (time[i] < tmin || time[i] > tmax) continue;
        *chi2 += resid[i];
        *ndata += 1;
    }

}

double kth_smallest(double a[], int n, int k) {
    // See: http://ndevilla.free.fr/median/median/
    int i,j,l,m ;
    double x, tmp ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                tmp = a[i];
                a[i] = a[j];
                a[j] = tmp;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

int single_fit (struct FitResult *fit, double *time, double *flux, 
    double *ferr, int npts, struct SearchResult *bump, 
    double outlier_sigma, double window) {
    
    double tstart,tstop,mean_flux,std_flux,tmax1,tmax2,max_flux;
    double *times_bump,*resid,pars[3],pars2[3],step[3],xmin[3],xmin2[3];
    double Fs,Fb,aux,chi2best,chi2best2,wyn0,wyn1,wyn2,wyn3,wyn4,wyn5;
    int i,nbad,nout,nc,nr,ne,*outlier,nbump,npt0,npt1,npt2,npt3,npt4,npt5;
    
    outlier = (int *) malloc(npts*sizeof(int));
    times_bump = (double *) malloc(npts*sizeof(double));
    resid = (double *) malloc(npts*sizeof(double));
    
    tstart = bump->ts-window;
    tstop = bump->tk+window;
    mean_flux = bump->mean;
    std_flux = bump->sigma;
    
    // Masking obvious outliers
    
    nbad = 0;
    for (i=0; i<npts; i++) {
        outlier[i] = 0;
        if ((time[i] < tstart || time[i] > tstop) && \
            (fabs(flux[i]-mean_flux) >= (outlier_sigma*std_flux))) {
            outlier[i] = 1;
            nbad += 1;
        }
    }
    
    // Finding initial t_0: maximal flux
    
    max_flux = 0.0;
    for (i=0; i<npts; i++) {
        if (outlier[i] == 1) continue;
        if (flux[i] > max_flux) {max_flux = flux[i]; tmax1 = time[i];}
    }
    
    // Finding initial t_0: median time of magnified data points
    
    nbump = 0;
    for (i=0; i<npts; i++) {
        if (time[i] < tstart || time[i] > tstop) continue;
        if (flux[i] >= (mean_flux+3.0*std_flux)) {
            times_bump[nbump] = time[i];
            nbump += 1;
        }
    }
    tmax2 = median(times_bump,nbump);
    
    // Perform an initial fit using tmax1 and tmax2 as starting points
    
    pars[0] = tmax1; pars[1] = 40.0; pars[2] = 0.5;
    pars2[0] = tmax2; pars2[1] = 40.0; pars2[2] = 0.5;
    step[0] = 10.0; step[1] = 20.0; step[2] = 0.5;
    
    nelmin(chi2,3,pars,xmin,&chi2best,0.01,step,50,10000,&nc,&nr,&ne,
        time,flux,ferr,outlier,npts);
    nelmin(chi2,3,pars2,xmin2,&chi2best2,0.01,step,50,10000,&nc,&nr,&ne,
        time,flux,ferr,outlier,npts);
    
    if (chi2best2 < chi2best) {
        xmin[0] = xmin2[0];
        xmin[1] = xmin2[1];
        xmin[2] = xmin2[2];
    }
    
    get_blending(xmin,time,flux,ferr,outlier,npts,&Fs,&Fb,resid,&aux);
    nout = remove_4_sigma_outliers(time,resid,outlier,npts);
    
    // Perform a 5-parameter fit
    
    nelmin(chi2,3,xmin,xmin,&chi2best,0.01,step,50,10000,&nc,&nr,&ne,
        time,flux,ferr,outlier,npts);
    
    xmin[1] = fabs(xmin[1]);
    xmin[2] = fabs(xmin[2]);
    
    get_blending(xmin,time,flux,ferr,outlier,npts,&Fs,&Fb,resid,&aux);
    
    calculate_chi2(&wyn0,&npt0,time,resid,outlier,npts,xmin[0]-xmin[1],
        xmin[0]+xmin[1]); // \chi^2 in the range |t-t0| < tE
    calculate_chi2(&wyn1,&npt1,time,resid,outlier,npts,xmin[0]-2*xmin[1],
        xmin[0]+2*xmin[1]); // \chi^2 in the range |t-t0| < 2*tE
    calculate_chi2(&wyn2,&npt2,time,resid,outlier,npts,xmin[0]-xmin[1],
        xmin[0]); // \chi^2 in the range t0-tE < t < t0
    calculate_chi2(&wyn3,&npt3,time,resid,outlier,npts,xmin[0],
        xmin[0]+xmin[1]); // \chi^2 in the range t0 < t < t0+tE
    calculate_chi2(&wyn4,&npt4,time,resid,outlier,npts,xmin[0]-2*xmin[1],
        xmin[0]); // \chi^2 in the range t0-2*tE < t < t0
    calculate_chi2(&wyn5,&npt5,time,resid,outlier,npts,xmin[0],
        xmin[0]+2*xmin[1]); // \chi^2 in the range t0 < t < t0+2*tE
    
    if (Fb >= BLENDING_LIMIT) {

        fit->n_outliers = nbad+nout;
        fit->t0 = xmin[0];
        fit->tE = xmin[1];
        fit->u0 = xmin[2];
        fit->Fs = Fs;
        fit->Fb = Fb;
        fit->error_code = ne;
        fit->chi2_all = chi2best;
        fit->n_all = npts-nout-nbad;
        fit->chi2_0 = wyn0;
        fit->n_0 = npt0;
        fit->chi2_1 = wyn2;
        fit->n_1 = npt2;
        fit->chi2_2 = wyn3;
        fit->n_2 = npt3;
        fit->chi2_3 = wyn1;
        fit->n_3 = npt1;
        fit->chi2_4 = wyn4;
        fit->n_4 = npt4;
        fit->chi2_5 = wyn5;
        fit->n_5 = npt5;

        free(outlier);
        free(times_bump);
        free(resid);
    
        return(0);
    }
        
    // If the 5-parameter fit yields unphysical results, perform
    // a 4-parameter fit
    
    nelmin(chi2_fixed,3,pars,xmin,&chi2best,0.01,step,50,10000,&nc,&nr,
        &ne,time,flux,ferr,outlier,npts);
    nelmin(chi2_fixed,3,pars2,xmin2,&chi2best2,0.01,step,50,10000,&nc,
        &nr,&ne,time,flux,ferr,outlier,npts);

    if (chi2best2 < chi2best) {
        xmin[0] = xmin2[0];
        xmin[1] = xmin2[1];
        xmin[2] = xmin2[2];
    }
        
    nelmin(chi2_fixed,3,xmin,xmin,&chi2best,0.01,step,50,10000,&nc,&nr,
        &ne,time,flux,ferr,outlier,npts);
    
    xmin[1] = fabs(xmin[1]);
    xmin[2] = fabs(xmin[2]);
    
    get_blending_fixed(xmin,time,flux,ferr,outlier,npts,&Fs,&Fb,resid,&aux);
    
    calculate_chi2(&wyn0,&npt0,time,resid,outlier,npts,xmin[0]-xmin[1],
        xmin[0]+xmin[1]); // \chi^2 in the range |t-t0| < tE
    calculate_chi2(&wyn1,&npt1,time,resid,outlier,npts,xmin[0]-2*xmin[1],
        xmin[0]+2*xmin[1]); // \chi^2 in the range |t-t0| < 2*tE
    calculate_chi2(&wyn2,&npt2,time,resid,outlier,npts,xmin[0]-xmin[1],
        xmin[0]); // \chi^2 in the range t0-tE < t < t0
    calculate_chi2(&wyn3,&npt3,time,resid,outlier,npts,xmin[0],
        xmin[0]+xmin[1]); // \chi^2 in the range t0 < t < t0+tE
    calculate_chi2(&wyn4,&npt4,time,resid,outlier,npts,xmin[0]-2*xmin[1],
        xmin[0]); // \chi^2 in the range t0-2*tE < t < t0
    calculate_chi2(&wyn5,&npt5,time,resid,outlier,npts,xmin[0],
        xmin[0]+2*xmin[1]); // \chi^2 in the range t0 < t < t0+2*tE
    
    fit->n_outliers = nbad+nout;
    fit->t0 = xmin[0];
    fit->tE = xmin[1];
    fit->u0 = xmin[2];
    fit->Fs = Fs;
    fit->Fb = Fb;
    fit->error_code = ne;
    fit->chi2_all = chi2best;
    fit->n_all = npts-nout-nbad;
    fit->chi2_0 = wyn0;
    fit->n_0 = npt0;
    fit->chi2_1 = wyn2;
    fit->n_1 = npt2;
    fit->chi2_2 = wyn3;
    fit->n_2 = npt3;
    fit->chi2_3 = wyn1;
    fit->n_3 = npt1;
    fit->chi2_4 = wyn4;
    fit->n_4 = npt4;
    fit->chi2_5 = wyn5;
    fit->n_5 = npt5;
    
    free(outlier);
    free(times_bump);
    free(resid);
    
    return(0);
    
}
