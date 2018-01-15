#include "util.h"
#include "phm.h"
#include "usage.h"

long readTable(gzFile f, double* X, double* bf, char* type, long* nexp, double** Xk, double** Bs, long nrow, long ncol, long* id1, long* id2, long* cumcol, long P, double sigma){
    // j all col in f
    // l only N & C
    // cumcol only N & C
    long i, j, k=0, l, m;
    long ldx = nrow*3 + P*2;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    double dcell;
    long lcell;
    double* bftemp; bftemp = (double*)calloc(10, sizeof(double));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        l=0;
        for(j=0; j<ncol; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t' || c=='\0'){
                    cell[k] = '\0';
                    if(type[j]=='B'){ //fprintf(stderr, "B");
                        if(nexp[j]==10){ //fprintf(stderr, "10\n");
                            sscanf(cell, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", bftemp, bftemp+1, bftemp+2, bftemp+3, bftemp+4, bftemp+5, bftemp+6, bftemp+7, bftemp+8, bftemp+9);
                        }else if(nexp[j]==6){ //fprintf(stderr, "6 %s\n", cell);
                            sscanf(cell, "%lf %lf %lf %lf %lf %lf", bftemp, bftemp+1, bftemp+2, bftemp+3, bftemp+4, bftemp+5);
                        }
                        for(m=0; m<nexp[j]; m++){
                            bf[i+nrow*m] = exp(bftemp[m]);
                            if(isnan(bf[i+nrow*m])>0){
                                bf[i+nrow*m] = 0.0;
                                //fprintf(stderr, "NaN is found at (%d, %d) element! (replaced by 0.0)\n", i, m);
                            }
                            //if(isinf(bf[i+nrow*m])>0){
                            //    bf[i+nrow*m] = 0.0;
                            //    fprintf(stderr, "Inf is found at (%d, %d) element! (replaced by 0.0)\n", i, m);
                            //}
                        }
                    }else if(type[j]=='O'){
                        sscanf(cell, "%lf", &dcell);
                        X[i*3+0+(2*P)*ldx] = X[i*3+1+(2*P)*ldx] = X[i*3+2+(2*P)*ldx] = dcell;
                    }else if(type[j]=='N'){
                        sscanf(cell, "%lf", &dcell);
                        X[(  cumcol[l])*ldx + i*3+0] = dcell;
                        X[(P+cumcol[l])*ldx + i*3+1] = dcell;
                        X[(P+cumcol[l])*ldx + i*3+2] = dcell;
                        for(m=0; m<nexp[j]; m++){
                            X[(  cumcol[l]+m+1)*ldx + i*3+0] = rk1(dcell, Xk[l][m]);
                            X[(P+cumcol[l]+m+1)*ldx + i*3+1] = rk1(dcell, Xk[l][m]);
                            X[(P+cumcol[l]+m+1)*ldx + i*3+2] = rk1(dcell, Xk[l][m]);
                        }
                        l++;
                    }else if(type[j]=='C'){
                        sscanf(cell, "%ld", &lcell);
                        if(lcell>nexp[j]-1) lcell = nexp[j]-1;
                        if(lcell>0){
                            X[i*3+0 + (  cumcol[l]+(lcell-1))*ldx] = 1.0;
                            X[i*3+1 + (P+cumcol[l]+(lcell-1))*ldx] = 1.0;
                            X[i*3+2 + (P+cumcol[l]+(lcell-1))*ldx] = 1.0;
                        }
                        l++;
                    }else if(type[j]=='J'){
                        sscanf(cell, "%ld", id1+i);
                    }else if(type[j]=='K'){
                        sscanf(cell, "%ld", id2+i);
                    }
                    k = 0;
                    break;
                }else{
                    cell[k++] = c;
                }
            }
        }
    }
    // penalty
    for(i=0; i<P*2; i++){// diag prior
        X[nrow*3+i+ldx*i] = 1.0/sigma;
    }
    // spline penalty
    l=0;
    for(j=0; j<ncol; j++){
        if(type[j]=='N'){
            for(i=0; i<nexp[j]; i++){
                for(m=0; m<nexp[j]; m++){
                    X[nrow*3+(  cumcol[l]+i+1)+(  cumcol[l]+m+1)*ldx] = Bs[l][i+m*nexp[j]];
                    X[nrow*3+(P+cumcol[l]+i+1)+(P+cumcol[l]+m+1)*ldx] = Bs[l][i+m*nexp[j]];
                }
            }
        }
        if(type[j]=='C' || type[j]=='N'){ l++; }
    }
    return 0;
}




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
    
    if(verbose>0) fprintf(stderr, "f1=%ld e1=%ld f2=%ld e2=%ld\n", f1[6], e1[6], f2[6], e2[6]);

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
                    //if(j==0){fprintf(stderr, "fwd pj=%lf qj=%lf %lf\n", p[j], q[j], Z[i+9*nrow]/pj);}
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
                    //if(j==0){fprintf(stderr, "rev pj=%lf qj=%lf %lf\n", p[j], q[j], Z[i+7*nrow]/pj);}
                    // finding most likely parent
                    if(Z[i+7*nrow] > maxparent[j]){ parent[j] = id1[i]; maxparent[j] = Z[i+7*nrow]; }
                }
            }
        }
    }
    return 0;
}



void RX(double* R, double* X, long H, long P, long ldx, double* Xt){
    long i, j, k;
    for(i=0; i<H; i++){// H : num of alt hypos
        for(j=0; j<P; j++){// P : num of covariates
            Xt[i+j*ldx] = 0.0;
            for(k=i; k<H; k++){
                Xt[i+j*ldx] += R[i+k*H] * X[k+j*ldx];
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
                            if(isnan(Y[i+(j-2)*nrow])>0){
                                Y[i+(j-2)*nrow] = 0.0;
                                fprintf(stderr, "NaN is found at (%d, %d) element! (replaced by 0.0)\n", i, j);
                            }
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
void MstepMultinom(double* Xt, double* X, double* R, double* z, double* beta, long N3, long P, long LDX){
    long i, j;
    long nsp, st;
    
    for(i=0; i<P; i++){
        for(j=0; j<P; j++){
            Xt[N3+i+j*LDX] = X[N3+i+j*LDX];
        }
    }
    
    // QR decomposition
    qr(Xt, R, LDX, P, LDX);
    // beta1 <- t(Xt) %*% y
    cblas_dgemv(CblasColMajor, CblasTrans, LDX, P, 1.0, Xt, LDX, z, 1, 0.0, beta+P+1, 1);
    // beta  <- backsolve(R, beta1)
    BackSolve(R, P, beta+P+1, beta);
}



double EstepMultinomSub(long a, long b, double* y, double* bf, double* X, double* beta, double* pphi, long N, long P, long LDX, long H, double* Xt, double* p, double* w, double* eta, double* z){
    int i, j, k;
    double lkhd=0.0;
    double wsq;
    // eta <- X %*% beta
    beta[P]=1.0;
    cblas_dgemv(CblasColMajor, CblasNoTrans, (b-a)*H, P+1, 1.0, X+a*H, LDX, beta, 1, 0.0, eta+a*H, 1);
    // p <- eta
    multisoftmax(eta+a*H, p+a*H, (b-a)*H, H);
    // y <- p, bf
    double tot;
    //double phiAll = 0.0;
    //double phi = 0.07862747919984960920381;// (*pphi); (*pphi) = 0.0;
    for(i=a; i<b; i++){
#ifdef FULLBF
        tot = (bf[i+0*N] + bf[i+1*N] + bf[i+2*N] + bf[i+3*N]) * (1.0 - p[i*H+0] - p[i*H+1] - p[i*H+2]);
        // phi is unknown
        //y[i*H+0] = (bf[i+4*N]*(1.-phi) + bf[i+5*N]*phi) * p[i*H+0];
        // phi was multiplied by bfs
        y[i*H+0] = (bf[i+4*N] + bf[i+5*N]) * p[i*H+0];
        
        y[i*H+1] = (bf[i+6*N] + bf[i+7*N]) * p[i*H+1];
        y[i*H+2] = (bf[i+8*N] + bf[i+9*N]) * p[i*H+2];
        tot += y[i*H+0] + y[i*H+1] + y[i*H+2];
        lkhd += log(tot);
        if(isnan(lkhd)>0){fprintf(stderr, "NaN produced at line %d\n", i+1); return 0.0;}
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
        for(j=0; j<H; j++){
            for(k=j; k<H; k++){
                z[i*H+j] += w[i*H*H+j+k*H]*(eta[i*H+k]-X[i*H+k+P*LDX]);
            }
        }
        RX(w+i*H*H, X+i*H, H, P, LDX, Xt+i*H);
    }
    free(ymp);
    return lkhd;
}

void* EstepMultinom(void* args){
    HIERARCHICAL2_MT* pmt = (HIERARCHICAL2_MT *)args;
    //fprintf(stderr ,"%ld %ld\n", pmt->a, pmt->b);
    pmt->lkhd = EstepMultinomSub(pmt->a, pmt->b, pmt->y, pmt->bf, pmt->X, pmt->beta, pmt->pphi, pmt->N, pmt->P, pmt->LDX, pmt->H, pmt->Xt, pmt->p, pmt->w, pmt->eta, pmt->z);
}


double vcovMultinom(double* X, double* Xt, double* R, double* y, double* p, long LDX, long N, long P, long H, double* res){
    double* V; V = (double*)calloc(P*2, sizeof(double));
    long i, j, k;
    for(i=0; i<N; i++){
        for(j=0; j<P; j++){
            Xt[i+j*N] = 0.0;
            for(k=0; k<H; k++){
                Xt[i+j*N] += X[i*H+k+j*LDX]*(y[i*H+k]-p[i*H+k]);
            }
        }
    }
    // QR decomposition
    qr(Xt, R, N, P, N);
    for(i=0; i<P; i++){
        clearAs0(V, 2*P);
        V[i]=1.0;
        ForwardSolve(R, P, V, V+P);
        BackSolve(R, P, V+P, V);
#ifdef PEN
        for(j=0; j<P; j++){ res[i*P+j] = V[j]; }
#else
        res[i] = sqrt(V[i]);
#endif
        //res[i]=sqrt(V[i]);
    }
    
    double vv = V[P-1];
    free(V);
    return vv;
    
}



// N = # all peak pairs
// P = # covariates
// X : N*3 x P
// Z : N x 10
void emATACMultinom(double* bf, double* X, long N, long P, long LDX, double* Z, double* z1, double* beta, long nthreads){
    long itr, i, j, k;
    
    //double* beta; beta = (double*)calloc(P*2, sizeof(double));
    double* eta;  eta  = (double*)calloc(N*3, sizeof(double));
    double* Xt;   Xt   = (double*)calloc(LDX*P, sizeof(double));
    double* w;    w    = (double*)calloc(N*3*3, sizeof(double));
    double* z;    z    = (double*)calloc(LDX, sizeof(double));// pseudo data
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
        pmt[tid].LDX  = LDX;
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
        
        MstepMultinom(Xt, X, R, z, beta, N*3, P, LDX);
        //fprintf(stderr, "phi = %lf\n beta=", phi);
        //fprintf(stderr, "beta=");printV2(beta, P);
        double totb=1.+exp(beta[0])+2.*exp(beta[1]);
        //fprintf(stderr, "pi=(%lf, %lf, %lf)\n", 1./totb, exp(beta[0])/totb, 2.*exp(beta[1])/totb);
        
    }
    // vcov
    vcovMultinom(X, Xt, R, y, p, LDX, N, P, 3L, beta+P);
    // Z
#ifdef FULLBF
    //double phi = 0.07862747919984960920381;
    double* p1; p1 = (double*)calloc(10, sizeof(double));
    for(i=0; i<N; i++){
        // independent
        p1[0] = p1[1] = p1[2] = p1[3] = (1.0 - p[i*3+0] - p[i*3+1] - p[i*3+2]);
        if(p1[0]<0.0){p1[0] = p1[1] = p1[2] = p1[3] = 0.0;}
        // phi is unknown
        //p1[4] = p[i*3+0]*(1.-phi);
        //p1[5] = p[i*3+0]*phi;
        // phi was multiplied by bf
        p1[4] = p[i*3+0];
        p1[5] = p[i*3+0];
        // causal
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


long parseCol(char* s, char** type, long** nexp){
    long i, n=1, m, mtot;
    for(i=0; i<strlen(s); i++){
        if(s[i]==','){n++;}
    }
    type[0] = (char*)calloc(n, sizeof(char));
    nexp[0] = (long*)calloc(n, sizeof(long));
    mtot = m = 0;
    for(i=0; i<n; i++){
        if(s[mtot+1]==','){
            type[0][i] = s[mtot]; m=1;
        }else{
            sscanf(s+mtot, "%c%ld%n", &(type[0][i]), &(nexp[0][i]), &m);
        }
        mtot += m + 1;
    }
    return n;
}

long countP(char** type, long** nexp, long n){
    long i, P=0;
    for(i=0; i<n; i++){
        if(type[0][i]=='N'){// numeric
            P += 1 + nexp[0][i];
        }else if(type[0][i]=='C'){// categorical
            P += nexp[0][i] - 1;
        }
    }
    return P;
}

int emColoc(double* bf, double* Z, long nrow, double* z1){
    long i, j, itr;
    double psi=0.1, del=0.08, lkhd=0.0, lkhd0=-1.0e20;
    double* d; d=(double*)calloc(nrow, sizeof(double));
    double* p; p=(double*)calloc(6, sizeof(double));
    for(itr=0; itr<1000; itr++){
        p[0]=p[1]=p[2]=p[3]=1.0-psi;
        p[4]=psi*(1.0-del);
        p[5]=psi*del;
        for(j=0; j<6; j++){ z1[j]=0.0; }
        for(i=0; i<nrow; i++){
            d[i]=0.0; // tot
            for(j=0; j<6; j++){
                d[i] += bf[i+j*nrow]*p[j];
            }
            if(d[i]>0.0 && isnan(d[i])==0){
                lkhd += log(d[i]);
                for(j=0; j<6; j++){
                    Z[i+j*nrow] = bf[i+j*nrow]*p[j]/d[i];
                    z1[j] += Z[i+j*nrow];
                    if(isinf(z1[j])>0){fprintf(stderr, "Inf at (%ld %ld)\n", i, j); return 1; }
                }
            }
        }
        psi = (z1[4]+z1[5])/(z1[0]+z1[1]+z1[2]+z1[3]+z1[4]+z1[5]);
        del = z1[5]/(z1[4]+z1[5]);
        if(isnan(psi)>0 || isnan(del)>0){printV2(z1, 6); fprintf(stderr, "Parameters become NaN!\n"); break;}
        fprintf(stderr, "[%ld] lkhd=%lf psi=%lf delta=%lf\n", itr, lkhd, psi, del);
        if(lkhd>lkhd0 && fabs(lkhd-lkhd0)<1.0e-6){
            fprintf(stderr, "Finished\n");
            break;
        }else{
            lkhd0 = lkhd;
            lkhd = 0.0;
        }
    }
    free(d);
    return 0;
}


int main(int argc, char** argv){
    long i, j, k, l;
    
    if(argc==1){usage_phm(); return 1;} 
    verbose = 0; for(i=0; i<argc; i++){if(strcmp(argv[i], "-v")==0){ verbose=1; }}
    
    long nthreads = 1;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--n-threads")==0 || strcmp(argv[k],"-t")==0){nthreads = (long)atoi(argv[k+1]);}}
    
    gzFile fi = NULL; // i = variant level
    gzFile fj = NULL; // j = feature level
    char* typi = NULL;
    long* nxpi = NULL;
    long* nxpj = NULL;
    long* cumcoli;
    long* id1;
    long* id2;
    long nvari;
    long nfeatures=0;
    long pi, Pi;// numbers of col and expanded col (p & P)
    long nrow=2;
    long ldx;
    
    long NB;
    
    double** Xki;
    double** Bsi;
    
    double* X;
    double* bf;
    
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-f")==0){ nfeatures = (long)atoi(argv[i+1]); }}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-r")==0){ nrow = (long)atoi(argv[i+1]); }}
    
    
    // input
    for(i=0; i<argc-1; i++){
        if(strcmp(argv[i], "-i")==0){
            fi = gzopen(argv[i+1], "rb6f");
            for(j=0; j<argc-1; j++){
                if(strcmp(argv[j], "-c")==0){
                    pi = parseCol(argv[j+1], &typi, &nxpi);
                    Pi = countP(&typi, &nxpi, pi)+1;
                    Xki = (double**)calloc(pi, sizeof(double*));
                    Bsi = (double**)calloc(pi, sizeof(double*));
                    cumcoli = (long*)calloc(Pi+1, sizeof(long));cumcoli[0]=1;
                    l=0;
                    for(k=0; k<pi; k++){
                        if(verbose>0)fprintf(stderr, "%ld %c %ld\n", Pi, typi[k], nxpi[k]);
                        if(typi[k]=='N' || typi[k]=='C'){// N or C
                            cumcoli[l+1] = cumcoli[l] + (typi[k]=='N' ? nxpi[k]+1 : nxpi[k]-1);
                            if(typi[k]=='N' && nxpi[k]>0){// only numeric and N of bspline bases > 0
                                getXk(&(Xki[l]),nxpi[k]);
                                getBS(&(Bsi[l]), Xki[l], nxpi[k]);
                                if(verbose>0){fprintf(stderr, "knots: "); printV2(Xki[l], nxpi[k]); }
                            }
                            if(verbose>0) fprintf(stderr, "cumcoli=%ld\n", cumcoli[l+1]);
                            l++;
                        }else if(typi[k]=='B'){NB = nxpi[k];}
                    }
                    break;
                }
            }
            nvari = l;
            if(typi==NULL){fprintf(stderr, "No column information for %s\n", argv[i+1]); return 1;}
            
            if(NB==10){
                ldx = nrow*3+Pi*2;
                X = (double*)calloc(ldx*(2*Pi+1), sizeof(double));
                bf= (double*)calloc(nrow*NB,      sizeof(double));
                id1 = (long*)calloc(nrow, sizeof(double));
                id2 = (long*)calloc(nrow, sizeof(double));
            
                for(j=0; j<nrow; j++){
                    X[0 *ldx + j*3+0] = 1.0;
                    X[Pi*ldx + j*3+1] = 1.0;
                    X[Pi*ldx + j*3+2] = 1.0;
                }
            
                readTable(fi, X, bf, typi, nxpi, Xki, Bsi, nrow, pi, id1, id2, cumcoli, Pi, 100.0);
                if(verbose>0){ 
                    //fprintf(stderr, "cumloci: "); printVL(cumloci, nfeatures+1);
                    //fprintf(stderr, "cumcoli: "); printVL(cumcoli, nvari+1);
                    //printM(X, nrow+Pi, Pi+1);
                    printM2(X+nrow*3, Pi*2, Pi*2+nrow*3, Pi*2);
                    //fprintf(stderr, "BF: "); printV(bf, 10);
                }
            }else if(NB==6){// Coloc
                ldx = nrow+Pi;
                //X = (double*)calloc(ldx*(Pi+1), sizeof(double));
                bf= (double*)calloc(nrow*NB,    sizeof(double));
                id1 = (long*)calloc(nrow, sizeof(double));
                id2 = (long*)calloc(nrow, sizeof(double));
                for(j=0; j<nrow; j++){
                    //X[0*ldx + j] = 1.0;
                }
                fprintf(stderr, "Loading data...");
                readTable(fi, NULL, bf, typi, nxpi, Xki, Bsi, nrow, pi, id1, id2, cumcoli, 0, 100.0);
                fprintf(stderr, "Done.\n\n");
            }
            break;
        }
    }
    if(fi==NULL){fprintf(stderr, "No input file!\n"); return 1;}
    
    double* Z; Z = (double*)calloc(nrow*NB, sizeof(double));
    double* z1; z1 = (double*)calloc(NB, sizeof(double));
    double* beta;
    if(NB==10){
        beta = (double*)calloc((Pi*2+1)*(Pi*2+1+1), sizeof(double));
        beta[0] =-4.007432;//-4.348672;
        beta[Pi]=-5.030170;//-4.166187;
        beta[Pi*2]=1.0;
        emATACMultinom(bf, X, nrow, Pi*2, ldx, Z, z1, beta, nthreads);
    }else if(NB==6){
        emColoc(bf, Z, nrow, z1);
    }
    
    
    //##########
    //# Output #
    //##########
    
    char* prefix;
    char* ofname;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-o")==0){ prefix = argv[i+1]; }}
    ofname = (char*)calloc(strlen(prefix)+100, sizeof(char));
    
    // directry exists?
    struct stat st = {0}; if (stat(prefix, &st) == -1) { mkdir(prefix, 0700); }
   
    if(NB==6){
        sprintf(ofname, "%s/pp.gz", prefix);
        gzFile outf = gzopen(ofname, "wb6f");
        for(i=0; i<nrow; i++){
            gzprintf(outf, "%ld\t%ld", id1[i], id2[i]);
            Z[i] += Z[i+4*nrow];
            for(j=0; j<4; j++){
                gzprintf(outf, "\t%lf", log(Z[i+j*nrow]));
            }
            gzprintf(outf, "\t%lf\n", log(Z[i+5*nrow]));
        }
        gzclose(outf);
        return 0;
    } 
    // hyper-parameters
    int binf=1;
    if(binf>0){// binary
        sprintf(ofname, "%s/peak_pair_level.bin", prefix);
        FILE* outf; outf = fopen(ofname, "wb");
#ifdef PEN
        fwrite(beta, sizeof(double), Pi*2+Pi*Pi*4, outf);
#else
        fwrite(beta, sizeof(double), Pi*4, outf);
#endif
        fclose(outf);
    }else{// gzipped
        sprintf(ofname, "%s/peak_pair_level.gz", prefix);
        gzFile outf = gzopen(ofname, "wb6f");
#ifdef PEN
        char sep = '\t';
        for(j=0; j<Pi*2+2; j++){
            sep='\t';
            for(i=0; i<Pi*2; i++){
                if(i==Pi*2-1){sep='\n';}
                gzprintf(outf, "%lf%c", beta[i+j*Pi*2], sep);
            }
        }
#else
        char sep = '\t';
        for(i=0; i<Pi*4; i++){
            if(i==Pi*4-1){sep='\n';}
            gzprintf(outf, "%lf%c", beta[i], sep);
        }
#endif
        gzclose(outf);
    }
    
    // posterior prob.
    binf=0;
    if(binf>0){// binary
        sprintf(ofname, "%s/pp.bin", prefix);
        FILE* outf; outf = fopen(ofname, "wb");
        fwrite(Z, sizeof(double), nrow*NB, outf);
        fclose(outf);
    }else{// gzipped
        sprintf(ofname, "%s/pp.gz", prefix);
        gzFile outf = gzopen(ofname, "wb6f");
        for(i=0; i<nrow; i++){
            gzprintf(outf, "%ld\t%ld", id1[i], id2[i]);
            Z[i] += Z[i+4*nrow] + Z[i+6*nrow] + Z[i+8*nrow];
            for(j=0; j<4; j++){
                gzprintf(outf, "\t%lf", log(Z[i+j*nrow]));
            }
            gzprintf(outf, "\t%lf", log(Z[i+5*nrow]));
            gzprintf(outf, "\t%lf", log(Z[i+7*nrow]));
            gzprintf(outf, "\t%lf\n", log(Z[i+9*nrow]));
        }
        gzclose(outf);
    }
    
    // PMR
    double* p;     p = (double*)calloc(nfeatures*2, sizeof(double)); // having 0 downstream peak
    double* q;     q = (double*)calloc(nfeatures*2, sizeof(double)); // having 0 upstream   peak
    double* tss; tss = (double*)calloc(nfeatures, sizeof(double)); // tss flag
    double* we;  we  = (double*)calloc(nfeatures, sizeof(double)); // we flag
    long*   parent;    parent    =   (long*)calloc(nfeatures, sizeof(long));   // parent id
    double* maxparent; maxparent = (double*)calloc(nfeatures, sizeof(double)); // maximum posterior prob for the "parent"
    
    //FILE* ftfbs; ftfbs = fopen("Data/segway.tss.pf.e.we.bin", "rb");
    //fseek(ftfbs, 8*nfeatures*0, SEEK_SET);
    //info = fread(tss, sizeof(double), nfeatures, ftfbs);
    //fseek(ftfbs, 8*nfeatures*3, SEEK_SET);
    //info = fread(we, sizeof(double), nfeatures, ftfbs);
    
    for(i=0; i<nfeatures; i++){
        tss[i]=1.0; 
        we[i]=1.0;
        parent[i] = -1;
    }
    getPeakwiseAnnot(id1, id2, Z, nrow, nfeatures, p, q, tss, we, parent, maxparent);
    
    sprintf(ofname, "%s/pmr.gz", prefix);
    gzFile outf_pmr = gzopen(ofname, "wb6f");
    for(i=0; i<nfeatures; i++){
        gzprintf(outf_pmr, "%lf\t%lf\t%lf\t%ld\t%lf\n", log((1.0-p[i])*(q[i])), p[i+nfeatures], q[i+nfeatures], parent[i]+1, log(maxparent[i]));
    }
    gzclose(outf_pmr);
    
    

    
}


int mainold(int argc, char** argv){
   
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
    
    
    long nrow=0, ncol=12, nfeatures;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-f")==0){ nfeatures = (long)atoi(argv[i+1]); }}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-r")==0){ nrow = (long)atoi(argv[i+1]); }}
    
    
    
    //long nrow=0, ncol=0;
    long binary_input = 0;
    /*for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--dimention")==0 || strcmp(argv[k],"-d")==0){
        sscanf(argv[k+1], "%ld,%ld", &nrow, &ncol);
        binary_input = 1;
        break;
    }}
    
    if(binary_input==0){// GZ files
        gzFile f = gzopen(argv[1], "rb6f");
        dim(f, &nrow, &ncol, 0);
        gzclose(f);
    }*/
    
    if(verbose>0) fprintf(stderr, "Data: %ld x %ld\n", nrow, ncol-2);
    double* X; X = (double*)calloc(nrow*(ncol-2), sizeof(double));
    double* Z; Z = (double*)calloc(nrow*(ncol-2), sizeof(double));
    double* z1; z1 = (double*)calloc(ncol-2, sizeof(double));
    long* id1; id1 = (long*)calloc(nrow, sizeof(long));
    long* id2; id2 = (long*)calloc(nrow, sizeof(long));
    
    
    
    
    
    // BF table read
    if(binary_input==0){// GZ
        gzFile f = NULL;
        for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-i")==0){f = gzopen(argv[i+1], "rb6f"); break;}}
        if(f==NULL){fprintf(stderr, "No input!\n"); return 1;}
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
   
    
    
    
    //##########
    //# Output #
    //##########
    
    char* prefix;
    char* ofname;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-o")==0){ prefix = argv[i+1]; }}
    ofname = (char*)calloc(strlen(prefix)+100, sizeof(char));
    
    // directry exists?
    struct stat st = {0}; if (stat(prefix, &st) == -1) { mkdir(prefix, 0700); }
    
    // hyper-parameters
    int binf=1;
    if(binf>0){// binary
        sprintf(ofname, "%s/peak_pair_level.bin", prefix);
        FILE* outf; outf = fopen(ofname, "wb");
#ifdef PEN
        fwrite(beta, sizeof(double), P+P*P, outf);
#else
        fwrite(beta, sizeof(double), P*2, outf);
#endif
        fclose(outf);
    }else{// gzipped
        sprintf(ofname, "%s/peak_pair_level.gz", prefix);
        gzFile outf = gzopen(ofname, "wb6f");
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
    }
    
    // posterior prob.
    binf=0;
    if(binf>0){// binary
        sprintf(ofname, "%s/pp.bin", prefix);
        FILE* outf; outf = fopen(ofname, "wb");
        fwrite(Z, sizeof(double), nrow*(ncol-2), outf);
        fclose(outf);
    }else{// gzipped
        sprintf(ofname, "%s/pp.gz", prefix);
        gzFile outf = gzopen(ofname, "wb6f");
        for(i=0; i<nrow; i++){
            gzprintf(outf, "%ld\t%ld", id1[i], id2[i]);
            Z[i] += Z[i+4*nrow] + Z[i+6*nrow] + Z[i+8*nrow];
            for(j=0; j<4; j++){
                gzprintf(outf, "\t%lf", log(Z[i+j*nrow]));
            }
            gzprintf(outf, "\t%lf", log(Z[i+5*nrow]));
            gzprintf(outf, "\t%lf", log(Z[i+7*nrow]));
            gzprintf(outf, "\t%lf\n", log(Z[i+9*nrow]));
        }
        gzclose(outf);
    }
    
    // PMR
    double* p;     p = (double*)calloc(nfeatures*2, sizeof(double)); // having 0 downstream peak
    double* q;     q = (double*)calloc(nfeatures*2, sizeof(double)); // having 0 upstream   peak
    double* tss; tss = (double*)calloc(nfeatures, sizeof(double)); // tss flag
    double* we;  we  = (double*)calloc(nfeatures, sizeof(double)); // we flag
    long*   parent;    parent    =   (long*)calloc(nfeatures, sizeof(long));   // parent id
    double* maxparent; maxparent = (double*)calloc(nfeatures, sizeof(double)); // maximum posterior prob for the "parent"
    
    //FILE* ftfbs; ftfbs = fopen("Data/segway.tss.pf.e.we.bin", "rb");
    //fseek(ftfbs, 8*nfeatures*0, SEEK_SET);
    //info = fread(tss, sizeof(double), nfeatures, ftfbs);
    //fseek(ftfbs, 8*nfeatures*3, SEEK_SET);
    //info = fread(we, sizeof(double), nfeatures, ftfbs);
    
    for(i=0; i<nfeatures; i++){
        tss[i]=1.0; 
        we[i]=1.0;
        parent[i] = -1;
    }
    getPeakwiseAnnot(id1, id2, Z, nrow, nfeatures, p, q, tss, we, parent, maxparent);
    
    sprintf(ofname, "%s/pmr.gz", prefix);
    gzFile outf_pmr = gzopen(ofname, "wb6f");
    for(i=0; i<nfeatures; i++){
        gzprintf(outf_pmr, "%lf\t%lf\t%lf\t%ld\t%lf\n", log((1.0-p[i])*(q[i])), p[i+nfeatures], q[i+nfeatures], parent[i]+1, log(maxparent[i]));
    }
    gzclose(outf_pmr);
        
    
 
}






