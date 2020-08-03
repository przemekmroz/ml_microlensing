struct SearchResult {
    double ts;
    double tk;
    double mean;
    double sigma;
    double chiout;
    double chibump;
    double chi2bump;
    int nbump;
    int ndet;
    double ampl;
    double duration;
};

int check_if_transient (double jd[], double mag[], double err[], int ngi, struct SearchResult *a, double tstart, double window, double step, double thresh_sigma, double outlier_sigma, int min_num_ampl_points);

void printSearchResult (struct SearchResult *a);
