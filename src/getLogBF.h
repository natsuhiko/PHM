
double getLogWABFMR2(double* g, double* x, int n, double* g2, double* y, int n2);
double getLogWABFMR(double* g, double* x, double* y, int n);
double getLogBF(double* g, double* y, int n, double sigma, double* work);
double getLogWABF(double* g, double* y, int n);
double getBeta(double* g, double* y, int n, double* se);

double getLogWABFInter(double* g, double* y, double* env, long N, double* work);

double getLogWABFfromBetaSE(double beta, double se);
