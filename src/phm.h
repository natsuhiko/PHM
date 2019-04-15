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

//#include <f2c.h>
//#include <blaswrap.h>
#include <clapack.h>
//#include <cblas.h>

#include <zlib.h>
#include <math.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


typedef struct{
        int tid;  // thread id
        //int npeaksperthread;  // number of peaks processed [tid*npeaks ... (tid+1)*npeaks-1]
        //int npeaks;  // total peaks
        //int nvars;
        //int* cumloci;
        //int* cumcats;
        //int P;
        //double* X;
        //double* BF;

        //double* beta;
        //double* Pi;

        //double* eta;
        //double* pjk;
        //double* z;
        //double* Z1;
        //double* w;
        //double* y;

        //double* Xt;
        //double* Ie;
        //double lkhd;
	
	int a;
	int b;
	double* y;
	double* bf;
	double* X;
	double* beta;
	double* pphi;
	int N;
	int P;
	int LDX;
	int H;
	double* Xt;
	double* p;
	double* w;
	double* eta;
	double* z;
	int offs;
	double lkhd;
}HIERARCHICAL2_MT;
