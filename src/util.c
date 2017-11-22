#include "util.h"


void printV(double* x, long n){
    long i;
    for(i=0; i<n; i++){
        fprintf(stdout, "%lf,", x[i]);
    }
    fprintf(stdout, "\n");
}

void printVL(long* x, long n){
    long i;
    for(i=0; i<n; i++){
        fprintf(stdout, "%ld,", x[i]);
    }
    fprintf(stdout, "\n");
}

void printM(double* x, long n, long m){
    long i, j;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            fprintf(stdout, "%lf,", x[i+j*n]);
        }
        fprintf(stdout, "\n");
    }
}

void printM2(double* x, long n, long ldx, long m){
    long i, j;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            fprintf(stderr, "%lf,", x[i+j*ldx]);
        }
        fprintf(stderr, "\n");
    }
}

int endWith(char* x, char* y){
    int i;
    int nx = strlen(x);
    int ny = strlen(y);
    for(i=0; i<ny; i++){
        if(x[nx-1-i]!=y[ny-1-i]){return 0;}
    }
    return 1;
}

void printV2(double* x, long n){
    long i;
    for(i=0; i<n; i++){
        fprintf(stderr, "%lf,", x[i]);
    }
    fprintf(stderr, "\n");
}


int mini(int a, int b){return a<b ? a : b;}
int maxi(int a, int b){return a>b ? a : b;}

double logit(double x){
    return log(x/(1.-x));
}
double expit(double x){
    return 1./(1.+exp(-x));
}
double rk1(double x, double z){
    return ((z-0.5)*(z-0.5)-1./12.)*((x-0.5)*(x-0.5)-1./12.)/4. - (pow(fabs(x-z)-0.5, 4.)-(fabs(x-z)-0.5)*(fabs(x-z)-0.5)/2. + 7./240.)/24.;
}
void rk(double* x0, double* x, double* xk, long n, long nk){
    long i, j;
    for(i=0; i<n; i++){
        if(x0[i]>=0.0){
            x[i] = x0[i];
            for(j=0; j<nk; j++){
                x[(j+1)*n+i] = rk1(x0[i], xk[j]);
            }
        }
    }
}



void expand(int* pos, double* val, int n, int* pos2, double* val2, int n2, double* w){// n < n2
    int i=0, j=0;
    while(i<n){
        if(pos2[j]<pos[i]){
            j++;
        }else if(pos[i]<pos2[j]){
            i++;
        }else if(pos[i]==pos2[j]){
            val2[j] = val[i];
            w[j] = 1.0;
            i++;
            j++;
        }
    }
}

void multisoftmax(double* eta, double* p, int N, int H){
    int i, j;
    double tot, geta;
    for(i=0; i<N/H; i++){
        geta = -100.0;
        for(j=0; j<H; j++){
            if(geta<eta[i*H+j]){geta=eta[i*H+j];}
        }
        geta -= 10.0;
        tot = exp(-geta);
        for(j=0; j<H; j++){
            tot += exp(eta[i*H+j]-geta);
        }
        for(j=0; j<H; j++){
            p[i*H+j] = exp(eta[i*H+j]-geta)/tot;
        }
    }
}

void softmax(double* eta, double* y, int n){
    int i;
    double offs=0.0, tot=0.0;
    for(i=0; i<n; i++){
        if(offs < eta[i]-20.0){offs = eta[i]-20.0;}
    }
    for(i=0; i<n; i++){
        y[i] = exp(eta[i] - offs);
        tot += y[i];
    }
    for(i=0; i<n; i++){
        y[i] /= tot;
    }
}

void wsoftmax(double* eta, double* y, double* w, int n){
    int i;
    double offs=0.0, tot=0.0;
    for(i=0; i<n; i++){if(w[i]>0.5){
        if(offs < eta[i]-20.0){offs = eta[i]-20.0;}
    }}
    for(i=0; i<n; i++){if(w[i]>0.5){
        y[i] = exp(eta[i] - offs);
        tot += y[i];
    }}
    for(i=0; i<n; i++){if(w[i]>0.5){
        y[i] /= tot;
    }else{y[i] = 0.0;}}
}

void clearAs0(double* x, int n){
    int i;
    for(i=0; i<n; i++){x[i]=0.0;}
}

void sumTo1(double* x, int n){
    double tot=0.0;
    int i;
    for(i=0; i<n; i++){tot += x[i];}
    for(i=0; i<n; i++){x[i] /= tot;}
}












long lmin(long a, long b){
    if(a>b){return b;}
    return a;
}



double nk_mean2(double* x, double* y, long n, long ldx){
    long i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx]*y[i*ldx] / dn;}
    return res;
}


double nk_mean(double* x, long n, long ldx){
    long i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx] / dn;}
    return res;
}

double nk_dsum(double* x, long n, long ldx){
    long i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

double nk_lsum2(double* x, double* p, long n, long ldx){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx]*log(p[i*ldx]);}
    return res;
}

double nk_lsum3(double* x, double* p, long n, long ldx){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx]*log(p[i*ldx]) + (1.-x[i*ldx])*log(1.-p[i*ldx]);}
    return res;
}

char dim(gzFile f, long* pN, long* pP, long skip){// N x P matrix
    
    long N=0;
    long P=0;
    long i;
    char sep;
    long ret;
    char c;
    while((c=gzgetc(f)) != EOF){
        if(N==0+skip && (c=='\n'||c==' '||c=='\t'||c==',')){
            if(c!='\n'){sep=c;} 
            P++;
        }
        if(c=='\n'){N++;}
    }
    gzseek(f, 0L, SEEK_SET);
    (*pN)=N;
    (*pP)=P;
    return sep;
}


long getXk(double** xk, long n){
    long i;
    xk[0] = (double*)calloc(n, sizeof(double));
    for(i=0; i<n; i++){
        xk[0][i] = ((double)i+1.-3./8.)/((double)n+1.-2.*3./8.);
    }
}

long getBS(double** bs, double* xk, long n){
    bs[0] = (double*)calloc(n*n, sizeof(double));
    long i, j;
    clearAs0(bs[0], n*n);
    for(i=0; i<n; i++){
        for(j=i; j<n; j++){
            bs[0][i+j*n] = rk1(xk[i], xk[j]);
        }
    }
    char uplo = 'U';
    long nn=n, lda=n, info;
    dpotrf_(&uplo, &nn, bs[0], &lda, &info);
    return info;
}


long setBS(double* xk, long n){
    BS = (double*)calloc(n*n, sizeof(double));
    long i, j;
    clearAs0(BS, n*n);
    for(i=0; i<n; i++){
        for(j=i; j<n; j++){
            BS[i+j*n] = rk1(xk[i], xk[j]);
        }
    }
    char uplo = 'U';
    long nn=n, lda=n, info;
    dpotrf_(&uplo, &nn, BS, &lda, &info);
    return info;
}
long setBS2(double* xk, long n){
    long np2=n+2;
    BS = (double*)calloc(np2*np2*np2*np2, sizeof(double));
    long i, j, k;
    clearAs0(BS, np2*np2*np2*np2);
    
    for(i=2; i<np2; i++){
        for(j=i; j<np2; j++){
            for(k=0; k<np2; k++){
                BS[i+j*np2*np2+k*np2*np2*np2+k*np2] += rk1(xk[i-2], xk[j-2]);
            }
        }
    }
    for(i=2; i<np2; i++){
        for(j=i; j<np2; j++){
            for(k=0; k<np2; k++){
                BS[i*np2+j*np2*np2*np2+k+k*np2*np2] += rk1(xk[i-2], xk[j-2]);
            }
        }
    }
    char uplo = 'U';
    long nn=n, lda=np2*np2, info;
    dpotrf_(&uplo, &nn, BS+(np2*np2*2+2), &lda, &info);
    nn=n, lda=np2*np2, info;
    dpotrf_(&uplo, &nn, BS+np2*np2*np2+np2+(np2*np2*2+2), &lda, &info);
    nn=np2*n, lda=np2*np2, info;
    dpotrf_(&uplo, &nn, BS+np2*np2*np2*2+np2*2, &lda, &info);
    
    for(i=0; i<np2*np2; i++){
        for(j=0; j<np2*np2; j++){
            if(i>j){ BS[i+j*np2*np2]=0.0; }
        }
    }
            
    return info;
}

void qr(double* A, double* R, long M, long N, long lda){
    long lwork = N*N;
    long info;
    double* work; work = (double*)calloc(lwork, sizeof(double));
    double* tau;  tau  = (double*)calloc(N,     sizeof(double));
    dgeqrf_(&M, &N, A, &lda, tau, work, &lwork, &info);
    //fprintf(stderr, "dgeqrf=%d\n", info);
    
    long K = N;
    long i,j;
    for(i=0; i<N; i++){
        for(j=i; j<N; j++){
            R[i+j*N] = A[i+j*lda];
        }
    }
    //fprintf(stderr, "%d %d %d %d %d %d\n", M, N, K, lda, lwork, info);
    dorgqr_(&M, &N, &K, A, &lda, tau, work, &lwork, &info);
    //fprintf(stderr, "dorgqr=%d\n", info);
    free(work);
    free(tau);
}


void ForwardSolve(double* R, long p, double* y, double* x){
    // To solve R^T x = y
    // ldy : leading dim of y
    long i, j;
    clearAs0(x, p);
    
    x[0] = y[0]/R[0];
    
    double tm = 0.0;
    for(i=1; i<p; i++){
        tm = 0.0;
        for(j=0; j<i; j++){
            tm += R[p*i+j]*x[j];
        }
        x[i] = (y[i]-tm)/R[p*i+i];
    }
}


void BackSolve(double* R, long p, double* y, double* x){
    // To solve y = R x
    long i, j;
    clearAs0(x, p);
    
    x[p-1] = y[p-1]/R[p*p-1];
    
    double tm = 0.0;
    for(i=p-2; i>=0; i--){
        tm = 0.0;
        for(j=p-1; j>i; j--){
            tm += R[p*j+i]*x[j];
        }
        x[i] = (y[i]-tm)/R[p*i+i];
    }
}





