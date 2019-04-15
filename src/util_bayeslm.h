#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

void coloc(double* lbfj, double* lbfk, double* eta0, double Pi1_j, double Pi1_k, double* w, int n, double* pp13);

void randomise(double* y, int n);
void randomise2(double* y, double* y2, int m, int n);

double rk1(double x, double z);
void rk(double* x0, double* x, double* xk, int n, int nk);

void printV(double* x, int n);
void printV2(double* x, int n);
void printVint(int* x, int n);
void printVlog(double* x, int n);

int bdfscanf1h(char* filename, double** x, int n0, int nskip);
int bdfscanf1(char* filename, double* x, int n0, int nskip);
int bdfscanfh(char* filename, double** x);
int gzfdscanf(char* filename, double** x);

double bdfsum(char* fname, int n, int skip);

double nk_summ(double* x, int n);
double nk_mean(double* x, int n);

double expit(double x);
void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double Pi1_a, double* w, int n, double* pp13, int* lci);

void expand(int* pos, double* val, int n, int* pos2, double* val2, int n2, double* w);
int mini(int a, int b);
int maxi(int a, int b);
void clearAs0(double* x, int n);
void softmax(double* eta, double* y, int n);
void wsoftmax(double* eta, double* y, double* w, int n);
//void pwhm(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp12, double* pphi);
void pwhm(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp12, double* pPhi0);
void pwhm12(double* bf1, double* bf2, double* bfmr1, double* bfrm2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss);
void pwhm12old(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss);
void pwhm1(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp5, double* pPhi0);
void pwhm13(double* lbf1, double* lbf2, double* lbfmr, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0);
void pwhm13old(double* lbf1, double* lbf2, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0);


void pwhmfm(double* lbfj, double* lbfk, double* lbfmrj, double* lbfmrk, double* etaa, double* etaj, double* etak, double Pi1_j, double Pi1_k, double* w, int n, double* Psi, double* Zjall);
void pwhmfm0(double* lbfj, double* etaj, double Pi1_j, double* w, int n, double* Zj0all);
//void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci);
void pwhmNewAtacGwas(double* lbf1, double* lbf2, double Pi1atac, double* w, int n, double* pp13);
void pwhmnewataceqtl(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j);

void pwhmnewataceqtlAllParam(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j);


