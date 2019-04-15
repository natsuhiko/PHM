#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
#include <cblas.h>

double* BS;
double* BT;

long verbose;

long endWith(char* x, char* y);
void printV(double* x, long n);
void printVL(long* x, long n);
void printM(double* x, long n, long m);
void printM2(double* x, long n, long ldx, long m);
void printV2(double* x, long n);
void printV2L(double* x, long n, long ldx);
void printV2n(double* x, long n);
void rk(double* x0, double* x, double* xk, long n, long nk);
void expand(long* pos, double* val, long n, long* pos2, double* val2, long n2, double* w);
long softmax(double* eta, double* y, long n);
void multisoftmax(double* eta, double* p, long N, long H);
void wsoftmax(double* eta, double* y, double* w, long n);
void clearAs0(double* x, long n);
void sumTo1(double* x, long n);
double logit(double x);
double expit(double x);
double rk1(double x, double z);
long mini(long a, long b);
long maxi(long a, long b);

long getXk(double** xk, long n);
long getBS(double** bs, double* xk, long n);
long setBS(double* xk, long n);
long setBS2(double* xk, long n);

long lmin(long a, long b);
double nk_mean(double* x, long n, long ldx);
double nk_mean2(double* x, double* y, long n, long ldx);
double nk_dsum(double* x, long n, long ldx);
double nk_lsum2(double* x, double* p, long n, long ldx);
double nk_lsum3(double* x, double* p, long n, long ldx);
char dim(gzFile f, long* pN, long* pP, long skip);

void qr(double* A, double* R, long M, long N, long lda);
void qr_lapack(double* A, double* R, long M, long N, long lda);
void ForwardSolve(double* R, long p, double* y, double* x);
void BackSolve(double* R, long p, double* y, double* x);

