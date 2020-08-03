struct FitResult {
    int n_outliers;
    double t0;
    double tE;
    double u0;
    double Fs;
    double Fb;
    int error_code;
    double chi2_all;
    int n_all;
    double chi2_0;
    int n_0;
    double chi2_1;
    int n_1;
    double chi2_2;
    int n_2;
    double chi2_3;
    int n_3;
    double chi2_4;
    int n_4;
    double chi2_5;
    int n_5;
};

double magnification (double pars[], double t);
void get_blending (double *params, double *time, double *flux, double 
    *f_err, int *outlier, int ndata, double *Fs, double *Fb, double 
    *resid, double *chi2);
void get_blending_fixed (double *params, double *time, double *flux, 
    double *f_err, int *outlier, int ndata, double *Fs, double *Fb, 
    double *resid, double *chi2);
double chi2 (double *params, double *time, double *flux, double *f_err, 
    int *outlier, int ndata);
double chi2_fixed (double *params, double *time, double *flux, double 
    *f_err, int *outlier, int ndata);
int remove_4_sigma_outliers (double *time, double *resid, int *flag, 
    int ndata);
void calculate_chi2 (double *chi2, int *ndata, double *time, double 
    *resid, int *outlier, int npts, double tmin, double tmax);
double kth_smallest(double a[], int n, int k);
int single_fit (struct FitResult *fit, double *time, double *flux, 
    double *ferr, int npts, struct SearchResult *bump, 
    double outlier_sigma, double window);
