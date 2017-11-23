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

#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
//#include <cblas.h>

#include <zlib.h>
#include <math.h>

typedef struct{
        long tid;  // thread id
        //long npeaksperthread;  // number of peaks processed [tid*npeaks ... (tid+1)*npeaks-1]
        //long npeaks;  // total peaks
        //long nvars;
        //long* cumloci;
        //long* cumcats;
        //long P;
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
	
	long a;
	long b;
	double* y;
	double* bf;
	double* X;
	double* beta;
	double* pphi;
	long N;
	long P;
	long LDXT;
	long H;
	double* Xt;
	double* p;
	double* w;
	double* eta;
	double* z;
	long offs;
	double lkhd;
}HIERARCHICAL2_MT;
