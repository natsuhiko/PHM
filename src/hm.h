#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <pthread.h>

//#include <cblas.h>
//#include <f2c.h>
//#include <clapack.h>

#include <zlib.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

double* beta0;
double* gamma0;

typedef struct{
	int tid;  // thread id
        int nfeaturesperthread;  // number of peaks processed [tid*npeaks ... (tid+1)*npeaks-1]
	int nfeatures;  // total peaks
	int* cumloci;
	int Pi;
	int Pj;
	double* X0;
	double* BF;

	double* beta;
	double* Pis;

	double* eta;
	double* pjk;
	double* z;
	double* Z1;
	double* w;
	double* y;
	
	double* Xt;
	double* Ie;
        double lkhd;
}HIERARCHICAL_MT;

