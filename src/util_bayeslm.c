#include "util_bayeslm.h"
#include "util.h"




// params : length of (2+P)N + (3P+2)P + 1
double lm(double* y, double* X, long N, long P, double* params){
    /*double* yhat;  yhat = (double*)calloc(N, sizeof(double));
    double* resid; resid = (double*)calloc(N, sizeof(double));
    double* Q;     Q = (double*)calloc(N*P, sizeof(double));
    double* R;     R = (double*)calloc(P*P, sizeof(double));
    double* H;     H = (double*)calloc(2*P*P, sizeof(double));*/
    double* beta;  beta = params;
    double* se;    se   = params+P;
    double* yhat;  yhat = params+2*P+1;
    double* resid; resid= params+2*P+1+N;
    double* Q;     Q    = params+2*P+1+N*2;
    double* R;     R    = params+2*P+1+N*2+N*P;
    double* H;     H    = params+2*P+1+N*2+N*P+P*P;
    double s2 = 0.0;
    long i;

    // Q R <- X
    nk_dcopy(X,Q,N*P);
    qr(Q, R, N, P, N);

    // beta <- R^-1 Q^T y
    nk_dgemvT(Q, N, P, N, y, H);
    BackSolve(R, P, H, beta);

    // yhat <- X %*% b, resid <- y - yhat, s2 = | resid |^2 / (N-P)
    nk_dgemv(X, N, P, N, beta, yhat);
    for(i=0; i<N; i++){ resid[i] = y[i] - yhat[i]; s2 += resid[i]*resid[i];}
    s2 /= (double)(N-P);

    // H <- R^-1 R^-T
    clearAs0(H, P*P*2);
    for(i=0; i<P; i++){
        H[P*i+i]=1.0;
        ForwardSolve(R, P, H+P*i, H+P*(P+i));
    }
    for(i=0; i<P; i++){
        BackSolve(R, P, H+P*(P+i), H+P*i);
        se[i] = sqrt(H[P*i+i]*s2);
    }
    params[P*2] = s2;
    /*free(Q);
    free(R);
    free(yhat);
    free(resid);
    free(H);*/
    return 0;
}


/*
double rk1(double x, double z){
    return ((z-0.5)*(z-0.5)-1./12.)*((x-0.5)*(x-0.5)-1./12.)/4. - (pow(fabs(x-z)-0.5, 4.)-(fabs(x-z)-0.5)*(fabs(x-z)-0.5)/2. + 7./240.)/24.;
}
void rk(double* x0, double* x, double* xk, int n, int nk){
    int i, j;
    for(i=0; i<n; i++){
        if(x0[i]>=0.0){
            x[i] = x0[i];
            for(j=0; j<nk; j++){
                x[(j+1)*n+i] = rk1(x0[i], xk[j]);
            }
        }
    }
}
*/

typedef struct{
    double val;
    int ord;
}ORDER;

int compare_ORDER( const void *c1, const void *c2 )
{
    ORDER test1 = *(ORDER *)c1;
    ORDER test2 = *(ORDER *)c2;
    
    double tmp1 = test1.val;   /* b を基準とする */
    double tmp2 = test2.val;
    
    if(tmp1-tmp2>0){return 1;}else if(tmp1-tmp2<0){return -1;}else{ return 0;}
}

void getRandomOrder(int n, int* ord){
    ORDER* array;
    int i;
    //srand((unsigned)(time(NULL)+getpid())); //srand( (unsigned int)time( NULL ) );
    array=(ORDER*)calloc(n, sizeof(ORDER));
    for(i=0; i<n; i++){
        array[i].val=rand();
        array[i].ord=ord[i];
        //fprintf(stderr, "%d %lf\n", array[i].ord, array[i].val);
    }
    qsort(array, n, sizeof(ORDER), compare_ORDER);
    for(i=0; i<n; i++){ord[i]=array[i].ord;}
    free(array);
}

void randomise(double* y, int n){
    double* y1;
    int* ord;
    int i;
    ord=(int*)calloc(n, sizeof(int));
    for(i=0; i<n; i++){ord[i]=i;}
    y1=(double*)calloc(n, sizeof(double));
    getRandomOrder(n, ord);
    for(i=0; i<n; i++){
        y1[i]=y[ord[i]];
    }
    for(i=0; i<n; i++){
        y[i]=y1[i];
    }
    free(y1);
    free(ord);
}

void randomise2(double* y, double* y2, int m, int n){
    double* y1;
    double* y21;
    int* ord;
    int i, j;
    ord=(int*)calloc(n, sizeof(int));
    for(i=0; i<n; i++){ord[i]=i;}
    y1=(double*)calloc(n, sizeof(double));
    y21=(double*)calloc(n*m, sizeof(double));
    getRandomOrder(n, ord);
    for(i=0; i<n; i++){
        y1[i]=y[ord[i]];
        for(j=0; j<m; j++){
            y21[i+j*n]=y2[ord[i]+j*n];
        }
    }
    for(i=0; i<n; i++){
        y[i]=y1[i];
        for(j=0; j<m; j++){
            y2[i+j*n]=y21[i+j*n];
        }
    }
    free(y1);
    free(ord);
}


/*
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

void printV(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stdout, "%lf ", x[i]);}
    fprintf(stdout, "%lf\n", x[i]);
}
void printV2(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf ", x[i]);}
    fprintf(stderr, "%lf\n", x[i]);
}
void printVint(int* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stdout, "%d ", x[i]);}
    fprintf(stdout, "%d\n", x[i]);
}
void printVlog(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf ", log(x[i]));}
    fprintf(stderr, "%lf\n", log(x[i]));
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


double nk_mean(double* x, int n){
    double res = 0.0;
    int i;
    for(i=0; i<n; i++){
        res += x[i];
    }
    return res / (double)n;
}
*/
double nk_sum(double* x, int n){
    double res = 0.0;
    int i;
    for(i=0; i<n; i++){
        res += x[i];
    }
    return res;
}


double bdfsum(char* fname, int n, int skip){
    //fprintf(stderr, "%d %d\n", n, skip);
    double res = 0.0;
    double* tmp; tmp = (double*)calloc(n, sizeof(double));
    int info, i;
    FILE* f;
    f = fopen(fname, "rb");
    fseek(f, skip*sizeof(double), SEEK_SET);
    info = fread(tmp, sizeof(double), n, f);
    fclose(f);
    res = nk_sum(tmp, n);
    //printV2(tmp,n);fprintf(stderr, "\n");
    free(tmp);
    return res;
}

/*
void sumTo1(double* x, int n){
    double tot=0.0;
    int i;
    for(i=0; i<n; i++){tot += x[i];}
    for(i=0; i<n; i++){x[i] /= tot;}
}
*/

int bdfscanf1h(char* filename, double** x, int n0, int nskip){
    FILE* f; f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    int n = ftell(f)/8 - nskip;
    fseek(f, nskip*sizeof(double), SEEK_SET);
    if(n0>0){n=n0;}
    (*x) = (double*)calloc(n, sizeof(double));
    n = fread(*x, sizeof(double), n, f);
    fclose(f);
    return n;
}

int bdfscanf1(char* filename, double* x, int n0, int nskip){
    FILE* f; f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    int n = ftell(f)/8 - nskip;
    fseek(f, nskip*sizeof(double), SEEK_SET);
    if(n0>0){n=n0;}
    n = fread(x, sizeof(double), n, f);
    fclose(f);
    return n;
}

int bdfscanfh(char* filename, double** x){
    return bdfscanf1h(filename, x, -1, 0);
}


int gzfdscanf(char* filename, double** x){// N x P matrix
    gzFile f = gzopen(filename, "rb6f");
    char c;
    int n=0;
    gzseek(f, 0L, SEEK_SET);
    int maxnc=0, nc=0;
    while((c=gzgetc(f)) != EOF){
        if(c=='\n'||c==' '||c=='\t'||c==','){
            n++;
            if(nc>maxnc){maxnc = nc;}
            nc=0;
        }else{
            nc++;
        }
    }
    //fprintf(stderr, "%d values read\n", n);
    int i=0;
    gzseek(f, 0L, SEEK_SET);
    (*x) = calloc(n, sizeof(double));
    char* cell; cell = (char*)calloc(maxnc+1, sizeof(char));
    nc=0;
    while((c=gzgetc(f)) != EOF){
        if(c=='\n'||c==' '||c=='\t'||c==','){
            cell[nc]='\0';
            sscanf(cell, "%lf", (*x)+i);
            i++;
            nc=0;
        }else{
            cell[nc++] = c;
        }
    }
    return n;
}


void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double Pi1_a, double* w, int n, double* pp13, int* lci){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    double Pi0_a = 1.0-Pi1_a;//-0.076026;
    //double Pi1_a = 1.0;//0.076026;
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;

    wsoftmax(eta0, p12,     w, n);
    wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_a       * exp(-2.*geta); 
    pp13[6] = Pi0_1       * exp(-2.*geta); 
    pp13[8] = Pi0_2       * exp(-2.*geta);
    for(j=0; j<n; j++){
      //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1_1 * Pi0_2 * p12[j+n]   * exp(lbf1[j]-geta); //
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); //
      //pp13[3] // linkage
        
      //pp13[4]  = Pi0_a;
        pp13[5] += Pi1_a * p12[j]     * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[6]  = Pi0_1;
        pp13[7] += Pi1_1 * p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta*2.); // 1 -> 2
      //pp13[8]  = Pi0_2;
        pp13[9] += Pi1_2 * p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta*2.); // 2 -> 1
        
        
        // same var
        diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j+n]*p12[j+n*2];
        /*if(lci[j]>=0){for(k=1; k<j+1; k++){// vars in the same peak
            if(lci[j]==lci[j-k]){
                diaglikel += Pi1_1*Pi1_2 *p12[j+n-k]*p12[j+n*2] * exp(lbf1[j-k]+lbf2[j] - geta*2.);
                diagprior += p12[j+n-k]*p12[j+n*2];
                diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2-k] * exp(lbf1[j]+lbf2[j-k] - geta*2.);
                diagprior += p12[j+n]*p12[j+n*2-k];
            }else{
                break;
            }
        }}*/
    }
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    if(pp13[3]<0.0){pp13[3]=0.0;}
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 10);
    
    free(p12);
}


// Param 5
void pwhmnewataceqtlAllParam(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    Pi1_1 = Pi1_2 = 1.0;
    double Pi0_1 = Pi1_1;
    double Pi0_2 = Pi1_2;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta0, p12,     w, n);
    //wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_2       * exp(-2.*geta);
    for(j=0; j<n; j++){
        pp13[1] += Pi1_1 * Pi0_2 * p12[j]     * exp(lbf1[j]-geta); // eqtl
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); // atac
        pp13[5] += Pi1_2 * p12[j+n*2]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio 
        diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];

    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}

void coloc(double* lbfj, double* lbfk, double* eta0, double Pi1_j, double Pi1_k, double* w, int n, double* pp13){
    int j, k;
    double* p12; p12 = (double*)calloc(n, sizeof(double));
    
    double Pi0_j = 1.0 - Pi1_j;
    double Pi0_k = 1.0 - Pi1_k;
    
    double Pi0_a=1.0, Pi1_a=1.0; 

    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbfj[j]){geta=lbfj[j];}
        if(geta<lbfk[j]){geta=lbfk[j];}
    }
    geta /= 2.0;
    geta -= 10.0;

    wsoftmax(eta0, p12, w, n);

    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_j*Pi0_k * exp(-2.*geta);
    pp13[4] = Pi0_a       * exp(-2.*geta);
    for(j=0; j<n; j++){
        pp13[1] += Pi1_j * Pi0_k * p12[j]          * exp(lbfj[j]-geta); // eqtl
        pp13[2] += Pi0_j * Pi1_k * p12[j]          * exp(lbfk[j]-geta); // atac
        pp13[5] += Pi1_a         * p12[j]          * exp(lbfj[j]+lbfk[j]  -geta*2.); // Pleio
        diaglikel += Pi1_j*Pi1_k * p12[j] * p12[j] * exp(lbfj[j]+lbfk[j] - geta*2.);
        diagprior += p12[j]*p12[j];

    }

    pp13[3] = (pp13[1]*pp13[2]/Pi0_j/Pi0_k - diaglikel)/(1.-diagprior); // linkage

    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);

    sumTo1(pp13, 6);

    free(p12);
}


void pwhmnewataceqtl(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta0, p12,     w, n);
    //wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_2       * exp(-2.*geta);
    //pp13[4] = (Pi0_1+Pi0_2)/2.0 * exp(-2.*geta);
    //pp13[4] = Pi0_1 * exp(-2.*geta);
    for(j=0; j<n; j++){
      //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1_1 * Pi0_2 * p12[j]     * exp(lbf1[j]-geta); // eqtl
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); // atac
      //pp13[3] // linkage
        
      //pp13[4]  = Pi0_0;
        pp13[5] += Pi1_2 * p12[j+n*2]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[5] += (Pi1_1+Pi1_2) * (p12[j]+p12[j+n*2])/4.0  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[5] += Pi1_1 * p12[j]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
        
        diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];
        
        // vars in the same peak
        /*if(lci[j]>=0){for(k=1; k<j+1; k++){
            //if(lci[j]==106){fprintf(stderr, "%d %d %s %d %s %d\n", j, j-k, rss[j], lci[j], rss[j-k], lci[j-k]);}
            if(lci[j]==lci[j-k]||lci[j-k]<0){
                if(fid2j==157756 &&strcmp("rs558245864",rss[j])==0 && strcmp("rs2409780",rss[j-k])==0){
                    fprintf(stderr, " %d %s %s %lf %lf %lf %lf\n", fid2j, rss[j-k],rss[j], p12[j-k], p12[j+n*2], p12[j], p12[j+n*2-k]);
                }
                //diaglikel += Pi1_1*Pi1_2 *p12[j-k]*p12[j+n*2] * exp(lbf1[j-k]+lbf2[j] - geta*2.);
                //diagprior += p12[j-k]*p12[j+n*2];
                //diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2-k] * exp(lbf1[j]+lbf2[j-k] - geta*2.);
                //diagprior += p12[j]*p12[j+n*2-k];
            }else{
                //if(strcmp("rs558245864",rss[j])==0||strcmp("rs558245864",rss[j-k])==0){fprintf(stderr, "  %s %s\n", rss[j-k],rss[j]);}
                break;
            }
        }}*/
    }
    
    //if(fid2j==157756)fprintf(stderr, "%lf %lf\n", diaglikel/(pp13[1]*pp13[2]/Pi0_1/Pi0_2), diagprior);
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}





void pwhmNewAtacGwas(double* lbf1, double* lbf2, double Pi1atac, double* w, int n, double* pp13){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* eta; eta = (double*)calloc(n,   sizeof(double));
    
    double Pi0atac = 1.0-Pi1atac;
    
    // Param 5
    Pi0atac = Pi1atac = 1.0;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta, p12,     w, n);
    wsoftmax(eta, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0atac * exp(-2.*geta);
    pp13[4] = Pi0atac * exp(-2.*geta);
    for(j=0; j<n; j++){
        //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1atac * p12[j]     * exp(lbf1[j]-geta); // atac
        pp13[2] += Pi0atac * p12[j+n*2] * exp(lbf2[j]-geta); // gwas
        //pp13[3] // linkage
        
        //pp13[4]  = Pi0_0;
        pp13[5] += Pi1atac * p12[j] * exp(lbf1[j]+lbf2[j] - geta*2.); // Pleio
        
        diaglikel += Pi1atac * p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];
    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0atac - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}







void pwhmfm0(double* lbfj, double* etaj, double Pi1_j, double* w, int n, double* Zj0all){
    int l;
    
    double* p12; p12 = (double*)calloc(n, sizeof(double));
    double* zj; zj = (double*)calloc(n, sizeof(double));
    
    double Pi0_j = 1.0-Pi1_j;
    
    
    double geta = 0;
    for(l=0; l<n; l++){
        if(geta<lbfj[l]){geta=lbfj[l];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(etaj, p12,     w, n);
    
    
    double totj = 0.0;
    // merginal for j
    for(l=0; l<n; l++){
        zj[l]  = Pi1_j * p12[l]   * exp(lbfj[l]-geta);
        totj  += zj[l];
    }
    
    for(l=0; l<n; l++){
        Zj0all[l] = zj[l];
    }
    
    free(p12);
    free(zj);
}




void pwhmfm(double* lbfj, double* lbfk, double* lbfmrj, double* lbfmrk, double* etaa, double* etaj, double* etak, double Pi1_j, double Pi1_k, double* w, int n, double* Psi, double* Zjall){
    int l;
    
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* zj; zj = (double*)calloc(n, sizeof(double));
    double* zk; zk = (double*)calloc(n, sizeof(double));
    double* zjk; zjk = (double*)calloc(n, sizeof(double));
    
    double Pi1_a = 0.07862747919984960920381;
    double Pi0_a = 1.0-Pi1_a;
    double Pi0_j = 1.0-Pi1_j;
    double Pi0_k = 1.0-Pi1_k;
    
    
    double geta = 0;
    for(l=0; l<n; l++){
        if(geta<lbfj[l]){geta=lbfj[l];}
        //if(geta<lbfk[l]){geta=lbfk[l];}
        //if(geta<lbfmrj[l]){geta=lbfmrj[l];}
        //if(geta<lbfmrk[l]){geta=lbfmrk[l];}
    }
    geta /= 2.0;
    //geta -= 10.0;
    
    wsoftmax(etaa, p12,     w, n);
    wsoftmax(etaj, p12+n,   w, n);
    wsoftmax(etak, p12+n*2, w, n);
    
    double offdiagprior=1.0;
    for(l=0; l<n; l++){
        offdiagprior -= p12[l+n]*p12[l+n*2];
    }
    
    double totj = 0.0;
    double totk = 0.0;
    // merginal for j and k
    for(l=0; l<n; l++){
        
        totj  += Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        totk  += Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        zj[l]  = Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        zk[l]  = Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        zjk[l]  = Psi[1] * Pi1_a * p12[l]     * exp(lbfj[l]+lbfk[l]  -geta*2.);
        zjk[l] += Psi[2] * Pi1_j * p12[l+n]   * exp(lbfj[l]+lbfmrj[l]-geta*2.);
        zjk[l] += Psi[3] * Pi1_k * p12[l+n*2] * exp(lbfk[l]+lbfmrk[l]-geta*2.);
        
    }
    
    double tmpzj;
    double tot = Psi[0] * Pi0_j*exp(-geta) * (Pi0_k*exp(-geta) + totk)      +      (Psi[1]*Pi0_a + Psi[2]*Pi0_j + Psi[3]*Pi0_k) * exp(-2.*geta);
    for(l=0; l<n; l++){
        tmpzj = zj[l];
        zj[l] = Psi[0] * zj[l]             * (Pi0_k*exp(-geta) + (totk - zk[l])/offdiagprior) + zjk[l];
        zk[l] = Psi[0] * zk[l]             * (Pi0_j*exp(-geta) + (totj - tmpzj)/offdiagprior) + zjk[l];
        
        tot += zj[l];
    }
    
    for(l=0; l<n; l++){
        Zjall[l] += zj[l];
        //zj[l] /= tot;
        //zk[l] /= tot;
    }
    
    free(p12);
    free(zjk);
    free(zj);
    free(zk);
}




















