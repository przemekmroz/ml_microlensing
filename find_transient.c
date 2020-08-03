#include <math.h>
#include <stdio.h>
#include "find_transient.h"

#define MAXOBS 20000
#define MIN_NUM_DATA 5

int check_if_transient (double jd[], double flux[], double ferr[], 
    int ngi, struct SearchResult *a, double tstart, double window, 
    double step, double thresh_sigma, double outlier_sigma, 
    int min_num_ampl_points) {

    int i,j,nwin,npts;
    double t1,t2,sum1,sum2,sum3,w;
    
    double mean,std,chiout;
    double thresh,chibump,chi2bump,dtbump,maxfluxbump,maxflux;
    int ndet,nbump,nd,nb;
    
    double t1_best=0.0, t2_best=0.0, av_best=0.0, std_best=0.0;
    double chi2bump_best=0.0, chiout_best=0.0, chibump_best=0.0;
    double dtbump_best=0.0, ampl=0.0, maxfluxbump_best=0.0;
    int nbump_best=0,ndet_best=0;
    
    nwin = (int)((jd[ngi-1]-tstart-window)/step)+1; // number of loops

    for (i=0; i<=nwin; i++) {
        t1 = tstart + (double)i*step;
        t2 = t1 + window;
        // Calculating mean brightness and scatter of data points
        // outside the window:
        sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; npts = 0;
        for (j=0; j<ngi; j++) {
            if (jd[j]<(t1-window) || jd[j]>(t2+window)) {
                w = 1.0/(ferr[j]*ferr[j]);
                sum1 += flux[j]*w;
                sum2 += flux[j]*flux[j]*w;
                sum3 += w;
                npts += 1;
            }
        }
        mean = sum1/sum3;
        std = sqrt(sum2*sum3-sum1*sum1)/sum3;
        chiout = std*std*sum3/((double)npts-1.0);
        // Removing 5 sigma outliers and recalculating mean flux
        // and scatter of data points:
        for (j=0; j<ngi; j++) {
            if (jd[j]<(t1-window) || jd[j]>(t2+window)) {
                if (fabs(flux[j]-mean) < outlier_sigma*std) continue;
                w = 1.0/(ferr[j]*ferr[j]);
                sum1 -= flux[j]*w;
                sum2 -= flux[j]*flux[j]*w;
                sum3 -= w;
                npts -= 1;
            }
        }
        //if (npts < MIN_NUM_DATA) continue; ///// pm
        mean = sum1/sum3;
        std = sqrt(sum2*sum3-sum1*sum1)/sum3;
        chiout = std*std*sum3/(double)npts;
        // Searching for consecutive points above the threshold inside 
        // the window
        thresh = mean + thresh_sigma*std;
        chi2bump = 0.0; chibump = 0.0; dtbump = 0.0; maxfluxbump = 0.0;
        nbump = 0; ndet = 0;
        sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; maxflux = 0.0;
        nb = 0; nd = 0;
        for (j=1; j<ngi; j++) {
            if (jd[j] < t1 || jd[j] > t2) continue;
            if (flux[j] >= thresh) {
                nb += 1;
                nd += 1; // detection on subtracted image
                sum1 += (flux[j]-mean)*(flux[j]-mean)/(ferr[j]*ferr[j]);
                sum2 += (flux[j]-mean)/ferr[j];
                sum3 += (jd[j]-jd[j-1]);
                if (flux[j] > maxflux) maxflux = flux[j];
            } else {
                if (sum1 > chi2bump && nb >= min_num_ampl_points) {
                    chi2bump = sum1; chibump = sum2; dtbump = sum3;
                    nbump = nb; maxfluxbump = maxflux; ndet = nd;
                }
                sum1 = 0.0; sum2 = 0.0; sum3 = 0.0;
                nb = 0; nd = 0; maxflux = 0.0;
            }
        }
        if (sum1 > chi2bump && nb >= min_num_ampl_points) {
            chi2bump = sum1; chibump = sum2; dtbump = sum3;
            nbump = nb; maxfluxbump = maxflux; ndet = nd;
        }
        /*
        printf("%10.5f %10.5f %6.3f %5.3f %7.2f %5d %7.1f %9.1f %4d\n",\
          t1, t2, mean, std, chiout, npts, chibump, chi2bump, nbump);
        */
        if (nbump >= min_num_ampl_points && chi2bump > chi2bump_best) {
            chi2bump_best = chi2bump;
            av_best = mean;
            std_best = std;
            chiout_best = chiout; 
            chibump_best = chibump;
            t1_best = t1;
            t2_best = t2;
            nbump_best = nbump;
            ndet_best = ndet;
            dtbump_best = dtbump;
            maxfluxbump_best = maxfluxbump;
        }
    }

    if (av_best > 0) ampl = 2.5*log10(maxfluxbump_best/av_best);
    else ampl = 0.0;

    a->ts = t1_best; a->tk = t2_best; a->mean = av_best; a->sigma = std_best;
    a->chiout = chiout_best; a->chibump = chibump_best; a->chi2bump = chi2bump_best;
    a->nbump = nbump_best; a->ndet = ndet_best; a->ampl = ampl; a->duration = dtbump_best;

    return(0);
}

void printSearchResult (struct SearchResult *a) {
    printf("%6.1f %6.1f %6.3f %5.3f %7.2f %7.1f %9.1f %4d %4d %5.2f %6.2f\n",\
    a->ts,a->tk,a->mean,a->sigma,a->chiout,a->chibump,a->chi2bump,a->nbump,a->ndet,a->ampl,a->duration);
}
