#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "find_transient.h"
#include "single.h"

#define MAXOBS 20000
#define MIN_PTS_NUMBER 30

/* Identifying microlensing events using deep neural nets
 * 
 * This program calculates the light curve features.
 * 
 * Usage:
 * ./get_features filename
 * 
 * filename should contain three columns: time, magnitude, and magnitude
 * error.
 * 
 * Compilation:
 * gcc -o get_features get_features.c find_transient.c single.c asa047.c -lm
 * 
 * P. Mroz @ Caltech, 1 Aug 2020
 */

int read_data_from_text_file (char *filename, double *time, double *mag, 
    double *magerr) {
    
    /* Reading data (time, magnitude, magnitude error) from a text file.
     * The input file should contain three columns. This function 
     * returns the number of data points in the file.
     */
    
    FILE *fp;
    int npts;
    
    npts = 0;
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    while (fscanf(fp,"%lf %lf %lf",&time[npts],&mag[npts],&magerr[npts])\
    != EOF) npts += 1;
    fclose(fp);
    
    return(npts);
}

void magnitudes_to_flux (double *time, double *flux, double *ferr, 
    double *mag, double *magerr, int npts) {

    /* Transforming magnitudes to flux */

    int i;
    
    for (i=0; i<npts; i++) {
        if (time[i] > 2450000.0) time[i] -= 2450000.0;
        flux[i] = pow(10.0,0.4*(18.0-mag[i]));
        ferr[i] = flux[i]*0.4*M_LN10*magerr[i];
    }
}

int delete_main_bump (double *time, double *flux, double *ferr, int npts, 
    struct SearchResult *bump) {

    int i,idx;
    double max_flux;
    
    // Finding event maximum:
    max_flux = 0.0;
    for (i=0; i<npts; i++) {
        if (time[i] < bump->ts || time[i] > bump->tk) continue;
        if (flux[i] > max_flux) { max_flux = flux[i]; idx = i; }
    }
    
    // Deleting the main bump:
    i = idx;
    while (i < npts && flux[i] > bump->mean) {
        flux[i] = -1.0;
        i += 1;
    }
    i += 1;
    while (i < npts && flux[i] > bump->mean) {
        flux[i] = -1.0;
        i += 1;
    }
    i = idx-1;
    while (i > 0 && flux[i] > bump->mean) {
        flux[i] = -1.0;
        i -= 1;
    }
    i -= 1;
    while (i > 0 && flux[i] > bump->mean) {
        flux[i] = -1.0;
        i -= 1;
    }
    
    // Cleanup data:
    
    idx = 0;
    for (i=0; i<npts; i++) {
        if (flux[i] < 0.0) continue;
        time[idx] = time[i];
        flux[idx] = flux[i];
        ferr[idx] = ferr[i];
        idx += 1;
    }
    
    return idx;
}

double clip (double x, double min, double max) {
    
    if (x < min) return min;
    else if (x > max) return max;
    return x;

}

void calculate_normalized_features (double *f, struct FitResult *fit, 
    struct SearchResult *bump, int flag) {
        
    f[0] = clip(log10(fit->tE+1.0e-5),-1.0,3.0)-1.3;
    f[1] = fit->u0;
    f[2] = fit->Fs/(fit->Fs+fit->Fb);
    f[3] = (double)fit->error_code-1.0;
    if (fit->n_all == 0) fit->n_all = 1;
    f[4] = clip(fit->chi2_all/fit->n_all,0.0,5.0)-1.0;
    if (fit->n_0 == 0) fit->n_0 = 1;
    f[5] = clip(fit->chi2_0/fit->n_0,0.0,5.0)-1.0;
    if (fit->n_1 == 0) fit->n_1 = 1;
    f[6] = clip(fit->chi2_1/fit->n_1,0.0,5.0)-1.0;
    if (fit->n_2 == 0) fit->n_2 = 1;
    f[7] = clip(fit->chi2_2/fit->n_2,0.0,5.0)-1.0;
    if (fit->n_3 == 0) fit->n_3 = 1;
    f[8] = clip(fit->chi2_3/fit->n_3,0.0,5.0)-1.0;
    if (fit->n_4 == 0) fit->n_4 = 1;
    f[9] = clip(fit->chi2_4/fit->n_4,0.0,5.0)-1.0;
    if (fit->n_5 == 0) fit->n_5 = 1;
    f[10] = clip(fit->chi2_5/fit->n_5,0.0,5.0)-1.0;
    f[11] = clip(18.0-2.5*log10(bump->mean),14.0,22.0)-18.0;
    f[12] = bump->chiout-1.0;
    f[13] = clip(bump->ampl,0.0,2.0)-1.0;
    f[14] = clip(log10(bump->duration),-0.5,2.5)-1.5;
    f[15] = (double)flag;
    
}

int main(int argc, char *argv[]) {
   
    double time[MAXOBS],mag[MAXOBS],magerr[MAXOBS],flux[MAXOBS],ferr[MAXOBS];
    double f[16];
    int res,flag,npts,npts2;
    
    double t_start = 8000.0;
    double window = 120.0;
    double step = 40.0;
    double thresh_sigma = 3.0;
    double outlier_sigma = 5.0;
    int min_num_ampl_points = 3;
    
    struct SearchResult bump,bump2;
    struct FitResult fit;
    
    FILE *fp;
    
    if (argc != 2) 
    {
        fprintf(stderr,"Usage: ./get_features filename\n");
        return(1);
    }

    /* Reading data from a text file */
    
    npts = read_data_from_text_file(argv[1],time,mag,magerr);
    if (npts <= MIN_PTS_NUMBER) return(0);

    /* Transforming magnitudes to flux */
    
    magnitudes_to_flux(time,flux,ferr,mag,magerr,npts);

    /* Check if transient */

    res = check_if_transient(time,flux,ferr,npts,&bump,t_start,window,
        step,thresh_sigma,outlier_sigma,min_num_ampl_points);

    if (res != 0 || bump.chibump < 32.0 || bump.chiout > 2.0) return(0);

    printf("%-20s %6.1f %6.1f %6.3f ",argv[1],bump.ts,bump.tk,bump.mean);
    printf("%5.3f %7.2f %7.1f ",bump.sigma,bump.chiout,bump.chibump);
    printf("%8.1f %4d %5.2f ",bump.chi2bump,bump.nbump,bump.ampl);
    printf("%5.2f ",bump.duration);
    
    /* PSPL model */
    
    res = single_fit(&fit,time,flux,ferr,npts,&bump,outlier_sigma,window);

    printf("%2d %9.4f %7.3f %6.4f ",fit.n_outliers,fit.t0,fit.tE,fit.u0);
    printf("%8.4f %8.4f %1d ",fit.Fs,fit.Fb,fit.error_code);
    printf("%7.2f %4d %7.2f %4d ",fit.chi2_all,fit.n_all,fit.chi2_0,fit.n_0);
    printf("%7.2f %4d %7.2f %4d ",fit.chi2_1,fit.n_1,fit.chi2_2,fit.n_2);
    printf("%7.2f %4d %7.2f %4d ",fit.chi2_3,fit.n_3,fit.chi2_4,fit.n_4);
    printf("%7.2f %4d ",fit.chi2_5,fit.n_5);

    /* Delete main bump */
    
    npts2 = delete_main_bump(time,flux,ferr,npts,&bump);
    
    /* Check if transient */
    
    res = check_if_transient(time,flux,ferr,npts2,&bump2,t_start,window,
          step,thresh_sigma,outlier_sigma,min_num_ampl_points);
    
    if (res == 0 && bump2.chibump > 32.0) flag = 1;
    else flag = 0;
        
    printf("%d\n",flag);
    
    /* Calculate normalized features and save to a text file */
    
    calculate_normalized_features(f,&fit,&bump,flag);
    
    fp = fopen("normalized_features.txt","a");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file normalized_features\n");
        return(0);
    }
    fprintf(fp,"%5.4f %5.4f %6.3f %4.1f ",f[0],f[1],f[2],f[3]);
    fprintf(fp,"%6.3f %6.3f %6.3f %6.3f ",f[4],f[5],f[6],f[7]);
    fprintf(fp,"%6.3f %6.3f %6.3f %6.3f ",f[8],f[9],f[10],f[11]);
    fprintf(fp,"%5.2f %5.2f %6.3f %3.1f\n",f[12],f[13],f[14],f[15]);
    fclose(fp);
    return(0);
}
