#include "util.h"
#include "phm.h"
#include "usage.h"

int getPeakwiseAnnot(long* id1, long* id2, double* Z, long nrow, long npeaks, double* p, double* q, double* tss, double* we, long* parent, double* maxparent){
    long pi, i, j, k;
    // 0 start
    for(i=0; i<nrow; i++){id1[i]--; id2[i]--;}
    // locations of peak j first seen
    long* f1; f1 = (long*)calloc(npeaks, sizeof(long));
    long* f2; f2 = (long*)calloc(npeaks, sizeof(long));
    long* e1; e1 = (long*)calloc(npeaks, sizeof(long));
    long* e2; e2 = (long*)calloc(npeaks, sizeof(long));
    long k1=0, k2=0;
    // init
    for(j=0; j<npeaks; j++){
        f1[j] = f2[j] = e1[j] = e2[j] = -1;
    }
    // finding start
    for(i=0; i<nrow; i++){
        if(f1[id1[i]]<0) f1[id1[i]] = i;
        if(f2[id2[i]]<0) f2[id2[i]] = i;
    }
    // finding end
    for(i=0; i<nrow; i++){
        e1[id1[i]] = i+1;
        e2[id2[i]] = i+1;
    }
    
    fprintf(stderr, "%ld %ld %ld %ld\n", f1[6], e1[6], f2[6], e2[6]);

    // calculating master/dependent prob
    double pj; // j is QTL
    for(j=0; j<npeaks; j++){
        p[j]=q[j] = 1.0;
        if(f1[j]>(-1)){// pair exists after j
            for(i=f1[j]; i<e1[j]; i++){
                if(id1[i]==j && we[id1[i]]>0.5 && tss[id2[i]]>0.5){
                    pj = Z[i+nrow] + Z[i+3*nrow] + Z[i+5*nrow] + Z[i+7*nrow] + Z[i+9*nrow];
                    p[j] *= (1.0-Z[i+7*nrow]/pj); // peak k is not a downstream peak for j | j is QTL
                    q[j] *= (1.0-(Z[i+5*nrow]+Z[i+9*nrow])/pj); // peak k is not a upstream   peak for j | j is QTL
                    p[j+npeaks] += Z[i+7*nrow]/pj;
                    q[j+npeaks] += Z[i+9*nrow]/pj;
                    if(j==0){fprintf(stderr, "fwd pj=%lf qj=%lf %lf\n", p[j], q[j], Z[i+9*nrow]/pj);}
                    // finding most likely parent
                    if(Z[i+9*nrow] > maxparent[j]){ parent[j] = id2[i]; maxparent[j] = Z[i+9*nrow]; }
                }
            }
        }
        if(f2[j]>(-1)){// pair exists before j
            for(i=f2[j]; i<e2[j]; i++){
                if(id2[i]==j && we[id2[i]]>0.5 && tss[id1[i]]>0.5){
                    pj = Z[i+2*nrow] + Z[i+3*nrow] + Z[i+5*nrow] + Z[i+7*nrow] + Z[i+9*nrow];
                    p[j] *= (1.0-Z[i+9*nrow]/pj); // peak k is not a downstream peak of j | j is QTL
                    q[j] *= (1.0-(Z[i+5*nrow]+Z[i+7*nrow])/pj); // peak k is not a upstream   peak of j | j is QTL
                    p[j+npeaks] += Z[i+9*nrow]/pj;
                    q[j+npeaks] += Z[i+7*nrow]/pj;
                    if(j==0){fprintf(stderr, "rev pj=%lf qj=%lf %lf\n", p[j], q[j], Z[i+7*nrow]/pj);}
                    // finding most likely parent
                    if(Z[i+7*nrow] > maxparent[j]){ parent[j] = id1[i]; maxparent[j] = Z[i+7*nrow]; }
                }
            }
        }
    }
    return 0;
}



void RX(double* R, double* X, long H, long P, long ldx, long ldxt, double* Xt){
    long i, j, k;
    for(i=0; i<H; i++){// H : num of alt hypos
        for(j=0; j<P; j++){// P : num of covariates
            Xt[i+j*ldxt] = 0.0;
            for(k=i; k<H; k++){
                Xt[i+j*ldxt] += R[i+k*H] * X[k+j*ldx];
            }
        }
    }
}


int makeW3(double* p, double* w){
    w[0] = sqrt(p[0]*(1.-p[0]));
    w[3] = -p[0]*p[1]/sqrt(p[0]*(1.-p[0]));
    w[4] = sqrt(p[1]*(1.-p[0]-p[1])/(1.-p[0]));
    w[6] = -p[0]*p[2]/sqrt(p[0]*(1.-p[0]));
    w[7] = -p[1]*p[2]/sqrt((1.-p[0])*p[1]*(1.-p[0]-p[1]));
    w[8] = sqrt(p[2]*(1.-p[0]-p[1]-p[2])/(1.-p[0]-p[1]));
    return 0;
}

long read2IDTable(gzFile f, long* id1, long* id2, double* Y, long nrow, long ncol){
    long i, j, k=0;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t' || c==' '){
                    cell[k] = '\0';
                    if(j==0){
                        sscanf(cell, "%ld", id1+i);
                    }else if(j==1){
                        sscanf(cell, "%ld", id2+i);
                    }else{
                        if(strcmp(cell,"-inf")==0){
                            Y[i+(j-2)*nrow] = 0.0;
                        }else{
                            sscanf(cell, "%lf", Y+i+(j-2)*nrow);
                            Y[i+(j-2)*nrow] = exp(Y[i+(j-2)*nrow]);
                        }
                    }
                    k = 0;
                    break;
                }else{
                    cell[k++] = c;
                }
            }
        }
    }
    return 0;
}

// Xt = W^1/2 %*% X
// z : pseudo data
void MstepMultinom(double* Xt, double* R, double* z, double* beta, long N3, long P, long LDXT, long offs){
    // penalty
    long i, j;
    double lambda = 1.0;
    long nsp, st;
    //// init
    
    //// fill with BS
#ifdef PEAKHEIGHT
    // init
    for(i=0; i<50; i++){for(j=0; j<P; j++){ Xt[N3+i+LDXT*j] = 0.0; }}
    // pleio
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 25, 15, 25, lambda, BS, 25, BT, 25, 0.0, Xt+N3, LDXT);
    // causal
    for(i=0; i<25; i++){
        for(j=0; j<25; j++){
            Xt[N3+(i+25)+LDXT*(j+15)] = lambda*BS[i + j*25];
        }
    }
    // print penalty mat
    //for(i=0; i<50; i++){for(j=0; j<P; j++){ fprintf(stderr, "%lf ", Xt[i+N3+j*LDXT]); } fprintf(stderr, "\n");}
#elif PEAKDIST
    // init 
    for(i=0; i<P; i++){for(j=0; j<P; j++){ Xt[N3+i+LDXT*j] = 0.0; }}
    for(nsp=0; nsp<2; nsp++){
        st=nsp*6+2;
        for(i=0; i<P; i++){
            for(j=0; j<P; j++){
                if(i>=st && j>=st && j>=i && i<st+4 && j<st+4){
                    Xt[N3+i+LDXT*j] = lambda*BS[i-st + (j-st)*4];
                }
            }
        }
    }
    //// print penalty mat
    //for(i=0; i<P; i++){for(j=0; j<P; j++){ fprintf(stderr, "%lf ", Xt[i+N3+j*LDXT]); } fprintf(stderr, "\n");}
#endif
    
    // QR decomposition
    qr(Xt, R, LDXT-offs, P-offs, LDXT);
    // beta1 <- t(Xt) %*% y
    cblas_dgemv(CblasColMajor, CblasTrans, LDXT-offs, P-offs, 1.0, Xt, LDXT, z, 1, 0.0, beta+P, 1);
    // beta  <- backsolve(R, beta1)
    BackSolve(R, P-offs, beta+P, beta);
}



double EstepMultinomSub(long a, long b, double* y, double* bf, double* X, double* beta, double* pphi, long N, long P, long LDXT, long H, double* Xt, double* p, double* w, double* eta, double* z, long offs){
    int i, j, k;
    double lkhd=0.0;
    double wsq;
    // eta <- X %*% beta
    beta[P]=1.0;
    cblas_dgemv(CblasColMajor, CblasNoTrans, (b-a)*H, P+1, 1.0, X+a*H, N*H, beta, 1, 0.0, eta+a*H, 1);
    // p <- eta
    multisoftmax(eta+a*H, p+a*H, (b-a)*H, H);
    // y <- p, bf
    double tot;
    //double phiAll = 0.0;
    double phi = 0.07862747919984960920381;// (*pphi); (*pphi) = 0.0;
    //phi = 0.05501029847975925229919;
    for(i=a; i<b; i++){
#ifdef FULLBF
        tot = (bf[i+0*N] + bf[i+1*N] + bf[i+2*N] + bf[i+3*N]) * (1.0 - p[i*H+0] - p[i*H+1] - p[i*H+2]);
        y[i*H+0] = (bf[i+4*N]*(1.-phi) + bf[i+5*N]*phi) * p[i*H+0];
        //y[i*H+0] = (bf[i+4*N] + bf[i+5*N]) * p[i*H+0];
        y[i*H+1] = (bf[i+6*N] + bf[i+7*N]) * p[i*H+1];
        y[i*H+2] = (bf[i+8*N] + bf[i+9*N]) * p[i*H+2];
        tot += y[i*H+0] + y[i*H+1] + y[i*H+2];
        lkhd += log(tot);
        y[i*H+0] /= tot;
        y[i*H+1] /= tot;
        y[i*H+2] /= tot;
#else
        tot = bf[i+0*N] * (1.0 - p[i*H+0] - p[i*H+1] - p[i*H+2]);
        y[i*H+0] = bf[i+1*N] * p[i*H+0];
        y[i*H+1] = bf[i+2*N] * p[i*H+1];
        y[i*H+2] = bf[i+3*N] * p[i*H+2];
        tot += y[i*H+0] + y[i*H+1] + y[i*H+2];
        lkhd += log(tot);
        y[i*H+0] /= tot;
        y[i*H+1] /= tot;
        y[i*H+2] /= tot;
#endif
    }
    
    // W <- p
    for(i=a; i<b; i++){ int info = makeW3(p+i*H, w+i*H*H); }
    
    // Xt <- W^1/2 %*% X
    // z  <- W^1/2 %*% X %*% beta + W^-1/2 %*% (y-p)
    double* ymp; ymp = (double*)calloc(H, sizeof(double));
    for(i=a; i<b; i++){
        for(j=0; j<H; j++){ ymp[j] = y[i*H+j] - p[i*H+j]; }
        ForwardSolve(w+i*H*H, H, ymp, z+i*H);
        if(offs==0){
            for(j=0; j<H; j++){
                for(k=j; k<H; k++){
                    z[i*H+j] += w[i*H*H+j+k*H]*(eta[i*H+k]-X[i*H+k+P*N*H]);
                }
            }
        }else{
            for(j=0; j<H; j++){
                for(k=j; k<H; k++){
                    z[i*H+j] += w[i*H*H+j+k*H]*(eta[i*H+k]-X[i*H+k+P*N*H]);
                }
            }
        }
        RX(w+i*H*H, X+i*H, H, P, N*H, LDXT, Xt+i*H);
    }
    free(ymp);
    return lkhd;
}

void* EstepMultinom(void* args){
    HIERARCHICAL2_MT* pmt = (HIERARCHICAL2_MT *)args;
    //fprintf(stderr ,"%ld %ld\n", pmt->a, pmt->b);
    pmt->lkhd = EstepMultinomSub(pmt->a, pmt->b, pmt->y, pmt->bf, pmt->X, pmt->beta, pmt->pphi, pmt->N, pmt->P, pmt->LDXT, pmt->H, pmt->Xt, pmt->p, pmt->w, pmt->eta, pmt->z, pmt->offs);
}


double vcovMultinom(double* X, double* Xt, double* R, double* y, double* p, long N, long P, long H, double* res, long offs){
    double* V; V = (double*)calloc(P*2, sizeof(double));
    long i, j, k;
    for(i=0; i<N/H; i++){
        for(j=0; j<P-offs; j++){
            Xt[i+j*(N/H)] = 0.0;
            for(k=0; k<H; k++){
                Xt[i+j*(N/H)] += X[i*H+k+j*N]*(y[i*H+k]-p[i*H+k]);
            }
        }
    }
    // QR decomposition
    qr(Xt, R, N/H, P-offs, N/H);
    for(i=0; i<P-offs; i++){
        clearAs0(V, 2*P);
        V[i]=1.0;
        ForwardSolve(R, P-offs, V, V+P);
        BackSolve(R, P-offs, V+P, V);
#ifdef PEN
        for(j=0; j<P-offs; j++){ res[i*(P-offs)+j] = V[j]; }
#else
        res[i] = sqrt(V[i]);
#endif
        //res[i]=sqrt(V[i]);
    }
    
    double vv = V[P-1-offs];
    free(V);
    return vv;
    
}



// N = # all peak pairs
// P = # covariates
// X : N*3 x P
// Z : N x 10
void emATACMultinom(double* bf, double* X, long N, long P, long LDXT, double* Z, double* z1, double* beta, long nthreads){
    long itr, i, j, k;
    
    //double* beta; beta = (double*)calloc(P*2, sizeof(double));
    double* eta;  eta  = (double*)calloc(N*3, sizeof(double));
    double* Xt;   Xt   = (double*)calloc(LDXT*P, sizeof(double));
    double* w;    w    = (double*)calloc(N*3*3, sizeof(double));
    double* z;    z    = (double*)calloc(LDXT, sizeof(double));// pseudo data
    double* p;    p    = (double*)calloc(N*3, sizeof(double));// pi
    double* y;    y    = (double*)calloc(N*3, sizeof(double));// pi
    double* R;    R    = (double*)calloc(P*P, sizeof(double));
    
    //double phi = 0.015649;
    double tot;
    double lkhd, lkhd0;
    lkhd0 = -1e10;
    long tid;
    long npp = N/nthreads + 1;
    if(verbose>0 && nthreads>1){fprintf(stderr, "Data split into %ld threads with %ld peak pairs each\n", nthreads, npp);}
    HIERARCHICAL2_MT* pmt; pmt=(HIERARCHICAL2_MT*)calloc(nthreads, sizeof(HIERARCHICAL2_MT));
    pthread_t* pid;       pid=(pthread_t*)calloc(nthreads, sizeof(pthread_t));
    for(tid=0; tid<nthreads; tid++){
        pmt[tid].tid  = tid;
        pmt[tid].a    = tid*npp;
        pmt[tid].b    = lmin((tid+1)*npp, N);
        pmt[tid].y    = y;
        pmt[tid].bf   = bf;
        pmt[tid].X    = X;
        pmt[tid].beta = beta;
        pmt[tid].N    = N;
        pmt[tid].P    = P;
        pmt[tid].LDXT = LDXT;
        pmt[tid].H    = 3;
        pmt[tid].Xt   = Xt;
        pmt[tid].p    = p;
        pmt[tid].w    = w;
        pmt[tid].eta  = eta;
        pmt[tid].z    = z;
    }
    for(itr=0; itr<10001; itr++){
        
        //lkhd = EstepATACMultinom(y, bf, X, beta, &phi, N, P, 3, Xt, p, w, eta, z, 0);
        for(tid=0; tid<nthreads; tid++){
            pmt[tid].lkhd = 0.0;
            //Estep(pmt);
            
            long pthflag;
            if( (pthflag = pthread_create(pid+tid, NULL, (void*)EstepMultinom, (void*)(pmt+tid))) != 0){
                fprintf(stderr, "Thread not created...aborted.\n");
                return;
            }
        }
        for(tid=0; tid<nthreads; tid++){
            long pthflag;
            if( (pthflag = pthread_join(pid[tid], NULL)) !=0 ){fprintf(stderr, "Thread not joined...aborted.\n"); return ;};
        }
        lkhd = 0.0;
        for(tid=0; tid<nthreads; tid++){lkhd += pmt[tid].lkhd;}
        
        
        
        
        
        if(isnan(lkhd)>0){fprintf(stderr, "NaN produced.\n"); return ;}
        if(fabs(lkhd-lkhd0)<1e-7){break;}else{lkhd0=lkhd; if(verbose>0){ fprintf(stderr, "[%ld] lkhd=%lf\n", itr, lkhd); } }
        
        MstepMultinom(Xt, R, z, beta, N*3, P, LDXT, 0);
        //fprintf(stderr, "phi = %lf\n beta=", phi);
        //fprintf(stderr, "beta=");printV2(beta, P);
        double totb=1.+exp(beta[0])+2.*exp(beta[1]);
        //fprintf(stderr, "pi=(%lf, %lf, %lf)\n", 1./totb, exp(beta[0])/totb, 2.*exp(beta[1])/totb);
        
    }
    // vcov
    vcovMultinom(X, Xt, R, y, p, N*3, P, 3L, beta+P, 0);
    // Z
#ifdef FULLBF
    double phi = 0.07862747919984960920381;
    //phi = 0.05501029847975925229919;
    double* p1; p1 = (double*)calloc(10, sizeof(double));
    for(i=0; i<N; i++){
        p1[0] = p1[1] = p1[2] = p1[3] = (1.0 - p[i*3+0] - p[i*3+1] - p[i*3+2]);
        if(p1[0]<0.0){p1[0] = p1[1] = p1[2] = p1[3] = 0.0;}
        p1[4] = p[i*3+0]*(1.-phi);
        p1[5] = p[i*3+0]*phi;
        p1[6] = p1[7] = p[i*3+1];
        p1[8] = p1[9] = p[i*3+2];
        tot = 0.0;
        for(j=0; j<10; j++){
            Z[i+j*N] = bf[i+j*N]*p1[j];
            tot += Z[i+j*N];
        }
        for(j=0; j<10; j++){
            Z[i+j*N] /= tot;
            z1[j] += Z[i+j*N];
        }
    }
    //printV(z1, 10);
#endif
    
}





int main(int argc, char** argv){
   
    long i, j, k;
    long info;
    
    if(argc==1){usage_phm(); return -1;}
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--verbose")==0 | strcmp(argv[k],"-v")==0){verbose=1;}}
    
    if(verbose>0){
        fprintf(stderr, "Command: ");
        for(i=0; i<argc; i++){ fprintf(stderr, "%s ", argv[i]); }
        fprintf(stderr, "\n\n");
    }
    
    long nthreads = 1;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--n-threads")==0 || strcmp(argv[k],"-t")==0){nthreads = (long)atoi(argv[k+1]);}}
    
    long nrow=0, ncol=0;
    long binary_input = 0;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--dimention")==0 || strcmp(argv[k],"-d")==0){
        sscanf(argv[k+1], "%ld,%ld", &nrow, &ncol);
        binary_input = 1;
        break;
    }}
    
    if(binary_input==0){// GZ files
        gzFile f = gzopen(argv[1], "rb6f");
        dim(f, &nrow, &ncol, 0);
        gzclose(f);
    }
    
    if(verbose>0) fprintf(stderr, "Data: %ld x %ld\n", nrow, ncol-2);
    double* X; X = (double*)calloc(nrow*(ncol-2), sizeof(double));
    double* Z; Z = (double*)calloc(nrow*(ncol-2), sizeof(double));
    double* z1; z1 = (double*)calloc(ncol-2, sizeof(double));
    long* id1; id1 = (long*)calloc(nrow, sizeof(long));
    long* id2; id2 = (long*)calloc(nrow, sizeof(long));
    
    
    // BF table read
    if(binary_input==0){// GZ
        gzFile f = gzopen(argv[1], "rb6f");
        read2IDTable(f, id1, id2, X, nrow, ncol);
    }else{// BIN
        FILE* f; f = fopen(argv[1], "rb");
        info = fread(X, sizeof(double), nrow*(ncol-2), f);
        fclose(f);
        f = fopen("Res2Splines/id1.bin", "rb");//hoge
        //f = fopen("id1.bin", "rb");
        info = fread(id1, sizeof(long), nrow, f);
        fclose(f);
        f = fopen("Res2Splines/id2.bin", "rb");
        //f = fopen("id2.bin", "rb");
        info = fread(id2, sizeof(long), nrow, f);
        fclose(f);
        //for(i=0; i<10; i++)fprintf(stderr, "%ld %ld\n", id1[i], id2[i]);
    }
    
    if(verbose>0) fprintf(stderr, "\nData was loaded...\n\n");
    
    // beta
    // coefs
#ifdef NIL
    long P = 2;
#elif PEAKHEIGHT
    long P = 40;
#elif PEAKDIST
    long P = 12;
#elif FIT
    long P = 5;
#endif
    
    // pairwise fit
    long fitpair = 0;
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--cooperative")==0){ fitpair=1; }}
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--collaborative")==0){ fitpair=2; }}
    //P += (fitpair==2 ? 1 : 0);
    
    
    
#ifdef PEN 
    double* beta; beta = (double*)calloc(P+P*P, sizeof(double));
#else
    double* beta; beta = (double*)calloc(2*P, sizeof(double));
#endif
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--input-beta")==0){
        FILE* fbeta = fopen(argv[k+1],"rb");
        info = fread(beta, sizeof(double), P, fbeta);
        fclose(fbeta);
        break;
    }}
    if(verbose>0){ 
        fprintf(stderr, "Init beta:"); 
        printV2(beta, P); 
        fprintf(stderr, "\n");
    }
    
    
   /* 
    // covariates
    FILE* fcovs; fcovs = fopen("Data/ph.bin", "rb");
    double* ph; ph = (double*)calloc(277128, sizeof(double));
    info = fread(ph, sizeof(double), 277128, fcovs);
    fclose(fcovs);
    fcovs = fopen("Data/midp.bin", "rb");
    double* midp; midp = (double*)calloc(277128, sizeof(double));
    info = fread(midp, sizeof(double), 277128, fcovs);
    fclose(fcovs);
    printV2(midp, 10);
    */
    
    
    
    
    
    ////creating cov mat
    double* Covs; Covs = (double*)calloc(nrow*3*(P+1), sizeof(double));
    
    
    
    
    
    
    
    
    // knots
    double apeakdis;
    double* xk6; xk6=(double*)calloc(4, sizeof(double)); xk6[0] = 0.1470588; xk6[1] = 0.3823529; xk6[2] = 0.6176471; xk6[3]=0.8529412;
    double* xk5; xk5=(double*)calloc(3, sizeof(double)); xk5[0] = 0.1923077; xk5[1] = 0.5000000; xk5[2] = 0.8076923;
    double* x1; x1 = (double*)calloc(5, sizeof(double));
    double* x2; x2 = (double*)calloc(5, sizeof(double));
    long jk2l[25] ={0,1, 2 ,3, 4, 
        1,5, 6, 7, 8, 
        2,6, 9,10,11, 
        3,7,10,12,13, 
        4,8,11,13,14};
    // penalty matrix for spline
#ifdef PEAKDIST
    setBS(xk6, 4);
#elif PEAKHEIGHT
    setBS2(xk5, 3);
    BT=(double*)calloc(15*25, sizeof(double));
    for(i=0; i<25; i++){BT[i+jk2l[i]*25]=1.0;}
    beta[40] = 1.0;
    FILE* covpeakdist; covpeakdist=fopen("Data/offs_peak_dist.bin","rb");
    info = fread(Covs+nrow*3*40, sizeof(double), nrow*3, covpeakdist);
    fclose(covpeakdist);
#elif FIT
    // tfbs
    
    long anrow=0, ancol=0;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--adim")==0){
        sscanf(argv[k+1], "%ld,%ld", &anrow, &ancol);
        break;
    }}
    if(anrow==0){fprintf(stderr, "No dimention info --adim\n"); return 1;}
    
    FILE* ftfbs=NULL;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--annot")==0){
        ftfbs = fopen(argv[k+1], "rb"); fprintf(stderr, "Annotation file: %s\n", argv[k+1]);
        break;
    }}
    if(ftfbs==NULL){fprintf(stderr, "No annotation --annot\n"); return 1;}
    
    long conp = 0;
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--conp")==0){ conp++; }}
    long symm = 1;
    long up=1;
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--downstream")==0){ symm=0; up=0;}}
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--upstream")==0){ symm=0; }}
    
    long tid=0;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--tid")==0){ tid =((long)atoi(argv[k+1]))-1; break;}}
    long tid2=tid/ancol; 
    tid = tid%ancol;
    fprintf(stderr, "tid=%ld tid2=%ld\n", tid+1, tid2+1);
    if(tid==tid2 && fitpair==2){P=4;}
    if(symm==0){P=3;}
    
    double* tfbs1;
    double* tfbs2;
    double avgt1, avgt2, avgt12;
    double* qhic; 
    if(anrow<nrow){
        tfbs1 = (double*)calloc(anrow, sizeof(double));
        tfbs2 = (double*)calloc(anrow, sizeof(double));
        fseek(ftfbs, 8*anrow*tid, SEEK_SET);
        info = fread(tfbs1, sizeof(double), anrow, ftfbs);
        fseek(ftfbs, 8*anrow*tid2, SEEK_SET);
        info = fread(tfbs2, sizeof(double), anrow, ftfbs);
        avgt1  = nk_mean(tfbs1, anrow, 1);
        avgt2  = nk_mean(tfbs2, anrow, 1);
        avgt12 = nk_mean2(tfbs1, tfbs2, anrow, 1);
        fprintf(stderr, "tfid1=%lf %lf %lf %lf %lf %lf...\n", tfbs1[0], tfbs1[1], tfbs1[2], tfbs1[3], tfbs1[4], tfbs1[5]);
        fprintf(stderr, "tfid2=%lf %lf %lf %lf %lf %lf...\n", tfbs2[0], tfbs2[1], tfbs2[2], tfbs2[3], tfbs2[4], tfbs2[5]);
    }else if(anrow==nrow){
        fprintf(stderr, "peak pair annot\n");
        qhic = (double*)calloc(nrow, sizeof(double));
        fseek(ftfbs, 8*nrow*tid, SEEK_SET);
        info = fread(qhic, sizeof(double), nrow, ftfbs);
    }
    
    
    // 0: no adj
    // 1: peak dist
    // 2: peak height
    // 3: peak dist & height
    long adj=0;
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--adj-pd")==0){ adj  = 1; break;}}
    for(k=0; k<argc; k++){if(strcmp(argv[k],"--adj-ph")==0){ adj += 2; break;}}
    
    fprintf(stderr, "Adj option: %ld\n", adj);
    
    
    
    //FILE* ftfbs; ftfbs = fopen("Data/hic_compartments6p2.bin", "rb"); fprintf(stderr, "HiC data\n"); long maxtf=8; // pair 8
    //FILE* ftfbs; ftfbs = fopen("Data/tfbs2.bin", "rb"); fprintf(stderr, "Encode ChIP data\n"); long maxtf=77; // pair 77
    //FILE* ftfbs; ftfbs = fopen("Data/segway.tss.pf.e.we.bin", "rb"); fprintf(stderr, "Segway data\n"); long maxtf=4; // pair 4
    //FILE* ftfbs; ftfbs = fopen("Data/histone2.bin", "rb"); fprintf(stderr, "Histone data\n"); long maxtf=11; // pair 11
    
    
    
    //FILE* fhicpeak; fhicpeak = fopen("Data/hicpeak.bin", "rb");
    //double* hicpeak; hicpeak = (double*)calloc(277128, sizeof(double));
    //info = fread(hicpeak, sizeof(double), 277128, fhicpeak);
    //FILE* fhicpeak = fopen("Data/hicpeakpair.bin", "rb");    // hic peak rao
    //FILE* fhicpeak = fopen("Data/hic_tad_intersect.bin", "rb");    // hic intersect rao
    //FILE* fhicpeak = fopen("Data/hic_tad_union.bin", "rb");    // hic union rao
    //FILE* fhicpeak = fopen("Data/chic.pr.ot.bin", "rb");    // chic pr pr
    //FILE* fhicpeak = fopen("Data/chic.pr.po.bin", "rb");    // chic pr other
    //qhic = 
    
    
    //betahat height
    double* betahat_ph; betahat_ph = (double*)calloc(40, sizeof(double));
    //betahat dist
    double* betahat;    betahat    = (double*)calloc(12, sizeof(double));
    FILE* fbetahat;
    if(adj==0){
        beta[0] = -4.3;// hoge
        beta[1] = -4.7;
        //beta[2] = 1.0;
    }else if(adj==1 || adj==3){
        fbetahat = fopen("Data/beta_peakdist.bin","rb");
        info = fread(betahat, sizeof(double), 12, fbetahat);
        fclose(fbetahat);
        printV2(betahat, 12);
    }
    if(adj==2){
        fbetahat = fopen("Data/beta40.bin","rb");
        info = fread(betahat_ph, sizeof(double), 40, fbetahat);
        fclose(fbetahat);
        printV2(betahat_ph, 40);
    }else if(adj==3){
        fbetahat = fopen("Data/beta40withPeakDistCorrected.bin","rb");
        info = fread(betahat_ph, sizeof(double), 40, fbetahat);
        fclose(fbetahat);
        printV2(betahat_ph, 40);
    }
    
    
#endif
    
    
    
    for(i=0; i<nrow; i++){
#ifdef NIL
        Covs[nrow*3*0 + i*3+0] = 1.0;
        Covs[nrow*3*1 + i*3+1] = 1.0;
        Covs[nrow*3*1 + i*3+2] = 1.0;
#elif PEAKHEIGHT
        x1[0] = x2[0] = 1.0;
        x1[1] = ph[id1[i]-1];
        x2[1] = ph[id2[i]-1];
        for(k=0; k<3; k++){
            x1[k+2] = rk1(x1[1], xk5[k]);
            x2[k+2] = rk1(x2[1], xk5[k]);
        }
        for(j=0; j<5; j++){
            for(k=0; k<5; k++){
                Covs[nrow*3*jk2l[j+k*5] + i*3+0] = (j==k ? x1[j]*x2[k] : x1[j]*x2[k]+x1[k]*x2[j]);
                Covs[nrow*3*(j+k*5+15)  + i*3+1] = x1[j]*x2[k];
                Covs[nrow*3*(j+k*5+15)  + i*3+2] = x1[k]*x2[j];
            }
        }
#elif PEAKDIST
        Covs[nrow*3*0 + i*3+0] = 1.0;
        Covs[nrow*3*6 + i*3+1] = 1.0;
        Covs[nrow*3*6 + i*3+2] = 1.0;
        
        apeakdis = fabs(midp[id1[i]-1]-midp[id2[i]-1])/500000.0;
        
        
        Covs[nrow*3*1 + i*3+0] = apeakdis;
        Covs[nrow*3*7 + i*3+1] = apeakdis;
        Covs[nrow*3*7 + i*3+2] = apeakdis;
        
        
        for(k=0; k<4; k++){
            Covs[nrow*3*(2+k) + i*3+0] = rk1(apeakdis, xk6[k]);
            Covs[nrow*3*(8+k) + i*3+1] = rk1(apeakdis, xk6[k]);
            Covs[nrow*3*(8+k) + i*3+2] = rk1(apeakdis, xk6[k]);
        }
#elif FIT
        // peak dist offset
        
        if(verbose>0){
           if(i==0)fprintf(stderr, "P=%ld, tid1=%ld, tid2=%ld, causal on pleio? %ld, fitpair? %ld\n", P, tid, tid2, conp, fitpair);
        }
        Covs[nrow*3*P + i*3+0] = betahat[0];
        Covs[nrow*3*P + i*3+1] = betahat[6];
        Covs[nrow*3*P + i*3+2] = betahat[6];
        
        apeakdis = fabs(midp[id1[i]-1]-midp[id2[i]-1])/500000.0;
        
        Covs[nrow*3*P + i*3+0] += apeakdis*betahat[1];
        Covs[nrow*3*P + i*3+1] += apeakdis*betahat[7];
        Covs[nrow*3*P + i*3+2] += apeakdis*betahat[7];
        
        for(k=0; k<4; k++){
            Covs[nrow*3*P + i*3+0] += rk1(apeakdis, xk6[k])*betahat[2+k];
            Covs[nrow*3*P + i*3+1] += rk1(apeakdis, xk6[k])*betahat[8+k];
            Covs[nrow*3*P + i*3+2] += rk1(apeakdis, xk6[k])*betahat[8+k];
        }
        
        // peak height offset
        
        x1[0] = x2[0] = 1.0;
        x1[1] = ph[id1[i]-1];
        x2[1] = ph[id2[i]-1];
        for(k=0; k<3; k++){
            x1[k+2] = rk1(x1[1], xk5[k]);
            x2[k+2] = rk1(x2[1], xk5[k]);
        }
        
        for(j=0; j<5; j++){
            for(k=0; k<5; k++){
                Covs[nrow*3*P + i*3+0] += (j==k ? x1[j]*x2[k] : (x1[j]*x2[k]+x1[k]*x2[j])/2.0)*betahat_ph[jk2l[j+k*5]];
                Covs[nrow*3*P + i*3+1] += x1[j]*x2[k] * betahat_ph[j+k*5+15];
                Covs[nrow*3*P + i*3+2] += x1[k]*x2[j] * betahat_ph[j+k*5+15];
            }
        }
        
        // intersept
        Covs[nrow*3*0 + i*3+0] = 1.0;
        Covs[nrow*3*1 + i*3+1] = 1.0;
        Covs[nrow*3*1 + i*3+2] = 1.0;
        
        if(anrow==nrow){
            // hic
            Covs[nrow*3*2 + i*3+0] = 0.0;
            Covs[nrow*3*2 + i*3+1] = qhic[i];
            Covs[nrow*3*2 + i*3+2] = qhic[i];
            Covs[nrow*3*3 + i*3+0] = qhic[i] + (conp>0 ? Covs[nrow*3*2 + i*3+0] : 0.0);
            Covs[nrow*3*3 + i*3+1] = 0.0     + (conp>0 ? Covs[nrow*3*2 + i*3+1] : 0.0);
            Covs[nrow*3*3 + i*3+2] = 0.0     + (conp>0 ? Covs[nrow*3*2 + i*3+2] : 0.0);
        }else{
            if(fitpair==0){
                // single
                Covs[nrow*3*2 + i*3+0] = 0.0;//2.0*avgt1;
                Covs[nrow*3*2 + i*3+1] = (up>0 ? tfbs1[id1[i]-1] : tfbs1[id2[i]-1]);
                Covs[nrow*3*2 + i*3+2] = (up>0 ? tfbs1[id2[i]-1] : tfbs1[id1[i]-1]);
                if(symm>0){
                    Covs[nrow*3*3 + i*3+0] = 0.0;//tfbs1[id1[i]-1] + tfbs1[id2[i]-1] + (conp>0 ? Covs[nrow*3*2 + i*3+0] : 0.0);
                    Covs[nrow*3*3 + i*3+1] = tfbs1[id2[i]-1]                   + (conp>0 ? Covs[nrow*3*2 + i*3+1] : 0.0);
                    Covs[nrow*3*3 + i*3+2] = tfbs1[id1[i]-1]                   + (conp>0 ? Covs[nrow*3*2 + i*3+2] : 0.0);
                }
            }else if(fitpair==1){
                // cooperative binding
                Covs[nrow*3*2 + i*3+0] = 0.0;//2.0*avgt12;
                Covs[nrow*3*2 + i*3+1] = (up>0 ? tfbs1[id1[i]-1]*tfbs2[id1[i]-1] : tfbs1[id2[i]-1]*tfbs2[id2[i]-1]);
                Covs[nrow*3*2 + i*3+2] = (up>0 ? tfbs1[id2[i]-1]*tfbs2[id2[i]-1] : tfbs1[id1[i]-1]*tfbs2[id1[i]-1]);
                if(symm>0){
                    Covs[nrow*3*3 + i*3+0] = 0.0;//tfbs1[id1[i]-1]*tfbs2[id1[i]-1] + tfbs1[id2[i]-1]*tfbs2[id2[i]-1] + (conp>0 ? Covs[nrow*3*2 + i*3+0] : 0.0);
                    Covs[nrow*3*3 + i*3+1] = tfbs1[id2[i]-1]*tfbs2[id2[i]-1]                                   + (conp>0 ? Covs[nrow*3*2 + i*3+1] : 0.0);
                    Covs[nrow*3*3 + i*3+2] = tfbs1[id1[i]-1]*tfbs2[id1[i]-1]                                   + (conp>0 ? Covs[nrow*3*2 + i*3+2] : 0.0);
                }
            }else{
                // collaborative binding
                Covs[nrow*3*2 + i*3+0] = tfbs1[id1[i]-1]*tfbs2[id2[i]-1] + tfbs2[id1[i]-1]*tfbs1[id2[i]-1];
                Covs[nrow*3*2 + i*3+1] = 0.0;//Covs[nrow*3*2 + i*3+0];
                Covs[nrow*3*2 + i*3+2] = 0.0;//Covs[nrow*3*2 + i*3+0];
                Covs[nrow*3*3 + i*3+0] = 0.0;
                Covs[nrow*3*3 + i*3+1] = tfbs1[id1[i]-1]*tfbs2[id2[i]-1];//Covs[nrow*3*2 + i*3+0];
                Covs[nrow*3*3 + i*3+2] = tfbs2[id1[i]-1]*tfbs1[id2[i]-1];//Covs[nrow*3*2 + i*3+0];
                if(tid!=tid2){
                    Covs[nrow*3*4 + i*3+0] = 0.0;
                    Covs[nrow*3*4 + i*3+1] = tfbs1[id1[i]-1]*tfbs2[id2[i]-1]*0.0 + tfbs2[id1[i]-1]*tfbs1[id2[i]-1];
                    Covs[nrow*3*4 + i*3+2] = tfbs1[id1[i]-1]*tfbs2[id2[i]-1]     + tfbs2[id1[i]-1]*tfbs1[id2[i]-1]*0.0;
                }
            }
        }
#endif

    }
    //cov output
    /*FILE* covpeakdist; covpeakdist=fopen("Data/offs_peak_height.bin","wb");
    fwrite(Covs+nrow*3*4, sizeof(double), nrow*3, covpeakdist);
    fclose(covpeakdist);*/
    if(verbose>0){
        for(i=0; i<3; i++){for(k=0; k<P+1; k++){fprintf(stderr, "%lf ", Covs[i+nrow*3*k]);} fprintf(stderr, "\n");} fprintf(stderr, "\n");
    }
    //fprintf(stderr, "means %lf %lf %lf %lf\n", nk_mean(Covs,nrow*3,1), nk_mean(Covs+nrow*3,nrow*3,1), nk_mean(Covs+nrow*6,nrow*3,1), nk_mean(Covs+nrow*9,nrow*3,1));
    
#ifdef NIL
    if(verbose>0) fprintf(stderr, "Null model...\n");
    beta[0]=-4.007432;//-4.348672;
    beta[1]=-5.030170;//-4.166187;
    emATACMultinom(X, Covs, nrow, P, nrow*3, Z, z1, beta, nthreads);
#elif PEAKHEIGHT
    fprintf(stderr, "Peak height...\n");
    emATACMultinom(X, Covs, nrow, P, nrow*3+50, Z, z1, beta, nthreads);
#elif PEAKDIST
    fprintf(stderr, "Peak dist...\n");
    emATACMultinom(X, Covs, nrow, P, nrow*3+P, Z, z1, beta, nthreads);
#elif FIT
    fprintf(stderr, "Fitted model...\n");
    emATACMultinom(X, Covs, nrow, P, nrow*3, Z, z1, beta, nthreads);
#endif
   
    
    
    
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--coefficients")==0 || strcmp(argv[k],"-c")==0){
        if(verbose>0)fprintf(stderr, "Output: %s\n", argv[k+1]);
        int binf=0;
        if(binf>0){
            FILE* outf; outf = fopen(argv[k+1], "wb");
#ifdef PEN
            fwrite(beta, sizeof(double), P+P*P, outf);
#else
            //printV(beta, P*2);
            fwrite(beta, sizeof(double), P*2, outf);
#endif
            fclose(outf);
            break;
        }else{
            gzFile outf = gzopen(argv[k+1], "wb6f");
#ifdef PEN
            char sep = '\t';
            for(j=0; j<P+1; j++){
                sep='\t';
                for(i=0; i<P; i++){
                    if(i==P-1){sep='\n';}
                    gzprintf(outf, "%lf%c", beta[i+j*P], sep);
                }
            }
#else
            char sep = '\t';
            for(i=0; i<P; i++){
                if(i==P-1){sep='\n';}
                gzprintf(outf, "%lf%c", beta[i], sep);
            }
#endif
            gzclose(outf);
            break;
        }
    }}
    
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--posterior-probability")==0 || strcmp(argv[k],"-p")==0){
        long kk, binf=0;
        for(kk=0; kk<argc; kk++){if(strcmp(argv[kk],"--bin")==0){binf++;}}
        if(binf>0){
            FILE* outf; outf = fopen(argv[k+1], "wb");
            fwrite(Z, sizeof(double), nrow*(ncol-2), outf);
            fclose(outf);
        }else{
            //for(i=0; i<10; i++)fprintf(stderr, "%ld %ld\n", id1[i], id2[i]);
            gzFile outf = gzopen(argv[k+1], "wb6f");
            for(i=0; i<nrow; i++){
                gzprintf(outf, "%ld\t%ld", id1[i], id2[i]);
                for(j=0; j<ncol-2; j++){
                    gzprintf(outf, "\t%lf", log(Z[i+j*nrow]));
                }
                gzprintf(outf, "\n");
            }
            gzclose(outf);
        }
        break;
    }}
    
#ifdef PEAKDIST
    long npeaks = 277128;
    double* p;     p = (double*)calloc(npeaks*2, sizeof(double)); // having 0 downstream peak
    double* q;     q = (double*)calloc(npeaks*2, sizeof(double)); // having 0 upstream   peak
    double* tss; tss = (double*)calloc(npeaks, sizeof(double)); // tss flag
    double* we;  we  = (double*)calloc(npeaks, sizeof(double)); // we flag
    long*   parent;    parent    =   (long*)calloc(npeaks, sizeof(long));
    double* maxparent; maxparent = (double*)calloc(npeaks, sizeof(double));
    
    FILE* ftfbs; ftfbs = fopen("Data/segway.tss.pf.e.we.bin", "rb");
    fseek(ftfbs, 8*npeaks*0, SEEK_SET);
    info = fread(tss, sizeof(double), npeaks, ftfbs);
    fseek(ftfbs, 8*npeaks*3, SEEK_SET);
    info = fread(we, sizeof(double), npeaks, ftfbs);
    
    for(i=0; i<npeaks; i++){
        tss[i]=1.0; 
        we[i]=1.0;
        parent[i] = -1;
    }
    fprintf(stderr, "tss=%lf we=%lf\n", tss[npeaks-1], we[npeaks-1]);
    
    getPeakwiseAnnot(id1, id2, Z, nrow, npeaks, p, q, tss, we, parent, maxparent);
    
    fprintf(stderr, "%lf %lf\n", p[0], q[0]);
    
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output-peak")==0){
        gzFile outf = gzopen(argv[k+1], "wb6f");
        for(i=0; i<npeaks; i++){
            //gzprintf(outf, "%lf\t%lf\t%lf\t%lf\t%ld\t%lf\n", log(1.0-p[i]), log(1.0-q[i]), log(p[i+npeaks]), log(q[i+npeaks]), parent[i]+1, log(maxparent[i]));
            gzprintf(outf, "%lf\t%lf\t%lf\t%lf\t%ld\t%lf\n", log(1.0-p[i]), log(q[i]), log(p[i+npeaks]), log(q[i+npeaks]), parent[i]+1, log(maxparent[i]));
        }
        gzclose(outf);
        break;
    }}
#endif
 
}






