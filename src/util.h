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

int verbose;

int endWith(char* x, char* y);
void printV(double* x, int n);
void printVL(int* x, int n);
void printM(double* x, int n, int m);
void printM2(double* x, int n, int ldx, int m);
void printV2(double* x, int n);
void printV2L(double* x, int n, int ldx);
void printV2n(double* x, int n);
void rk(double* x0, double* x, double* xk, int n, int nk);
void expand(int* pos, double* val, int n, int* pos2, double* val2, int n2, double* w);
int softmax(double* eta, double* y, int n);
void multisoftmax(double* eta, double* p, int N, int H);
void wsoftmax(double* eta, double* y, double* w, int n);
void clearAs0(double* x, int n);
void sumTo1(double* x, int n);
double logit(double x);
double expit(double x);
double rk1(double x, double z);
int mini(int a, int b);
int maxi(int a, int b);

int getXk(double** xk, int n);
int getBS(double** bs, double* xk, int n);
int setBS(double* xk, int n);
int setBS2(double* xk, int n);

int lmin(int a, int b);
double nk_mean(double* x, int n, int ldx);
double nk_mean2(double* x, double* y, int n, int ldx);
double nk_dsum(double* x, int n, int ldx);
double nk_lsum2(double* x, double* p, int n, int ldx);
double nk_lsum3(double* x, double* p, int n, int ldx);
char dim(gzFile f, int* pN, int* pP, int skip);

void qr(double* A, double* R, int M, int N, int lda);
void qr_lapack(double* A, double* R, int M, int N, int lda);
void ForwardSolve(double* R, int p, double* y, double* x);
void BackSolve(double* R, int p, double* y, double* x);

