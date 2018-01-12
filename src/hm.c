#include "hm.h"
#include "util.h"
#include "usage.h"

int cmpdr (const void * a, const void * b){
    if( *(double*)b - *(double*)a > 0.0){
        return 1;
    }else if(*(double*)b - *(double*)a < 0.0){
        return -1;
    }else{
        return 0;
    }
}

// P = expanded ncol
long readTable(gzFile f, double* X, double* bf, char* type, long* nexp, double** Xk, double** Bs, long nrow, long ncol, long* cumloci, long* cumcol, long P, double sigma){
    // j all col in f
    // l only N & C
    // cumcol only N & C
    long i, j, k=0, l, m;
    long ldx = nrow + P;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    double dcell;
    long lcell;
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        l=0;
        for(j=0; j<ncol; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t' || c=='\0'){
                    cell[k] = '\0';
                    if(type[j]=='B'){
                        sscanf(cell, "%lf", bf+i);
                        bf[i] = exp(bf[i]);
                        if(isnan(bf[i])>0){fprintf(stderr, "Bayes factor is NaN at %ld row.\n", i+1);}
                        if(isinf(bf[i])>0){fprintf(stderr, "Bayes factor is Inf at %ld row.\n", i+1); bf[i]=exp(500.0);}
                    }else if(type[j]=='O'){
                        sscanf(cell, "%lf", X+i+P*ldx);
                    }else if(type[j]=='N'){
                        sscanf(cell, "%lf", &dcell);
                        X[i+cumcol[l]*ldx] = dcell;
                        for(m=0; m<nexp[j]; m++){
                            X[i+(cumcol[l]+m+1)*ldx] = rk1(dcell, Xk[l][m]);
                        }
                        l++;
                    }else if(type[j]=='C'){
                        sscanf(cell, "%ld", &lcell);
                        if(lcell>nexp[j]-1) lcell = nexp[j]-1;
                        if(lcell>0) X[i+(cumcol[l]+(lcell-1))*ldx] = 1.0;
                        l++;
                    }else if(type[j]=='I'){
                        sscanf(cell, "%ld", &lcell);
                        cumloci[lcell]=i+1;
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
    for(i=0; i<P; i++){// diag prior
        X[nrow+i+ldx*i] = 1.0/sigma;
    }
    // spline penalty
    l=0;
    for(j=0; j<ncol; j++){
        if(type[j]=='N'){
            for(i=0; i<nexp[j]; i++){
                for(m=0; m<nexp[j]; m++){
                    X[nrow+(cumcol[l]+i+1)+(cumcol[l]+m+1)*ldx] = Bs[l][i+m*nexp[j]];
                }
            }
        }
        if(type[j]=='C' || type[j]=='N'){ l++; }
    }
    return 0;
}


void dgemv(double* A, long N, long M, long lda, double* x, double* y){
    long i, j;
    for(i=0; i<N; i++){
        y[i] = 0.0;
        for(j=0; j<M; j++){
            y[i] += A[i+j*lda]*x[j];
        }
    }
}


void dgemvT(double* A, long N, long M, long lda, double* x, double* y){
    long i, j;
    for(j=0; j<M; j++){
        y[j] = 0.0;
        for(i=0; i<N; i++){
            y[j] += A[i+j*lda]*x[i];
            //if(isnan(y[j])>0){fprintf(stderr, "NaN generated in dgevmT at (%ld, %ld)\n", i, j);}
        }
    }
}


void* Estep(void *args){
    HIERARCHICAL_MT* pmt = (HIERARCHICAL_MT *)args;
    
    long tid = pmt->tid;
    long nfeaturesperthread = pmt->nfeaturesperthread;
    
    long P = pmt->Pi;
    
    long nfeatures = pmt->nfeatures;
    long* cumloci; cumloci = pmt->cumloci;
    long L = cumloci[nfeatures];
    
    long ldx = L+P;
    
    //fprintf(stderr, "L=%d P=%d nvars=%d tid=%d nppt=%d\n", L, P, nvars, tid, nfeaturesperthread);
    
    double* X0; X0 = pmt->X0;
    double* BF; BF = pmt->BF;
    
    double* beta; beta = pmt->beta;
    double* Pi;   Pi   = pmt->Pis;
    
    double* eta; eta = pmt->eta;
    double* pjk; pjk = pmt->pjk;
    double* z;   z   = pmt->z;
    double* Z1;  Z1  = pmt->Z1;
    double* w;   w   = pmt->w;
    double* y;   y   = pmt->y;
    double* Xt;  Xt  = pmt->Xt;
    double* Ie;  Ie  = pmt->Ie;  clearAs0(Ie, P*P);
    
    double tot;
    
    double* xb; xb = (double*)calloc(P, sizeof(double));
    double* XjTZj; XjTZj = (double*)calloc(P, sizeof(double));
    
    long j, k, l, l2, pstart = tid*nfeaturesperthread, pend, nanflag=0;;
    pend = ((tid+1)*nfeaturesperthread < nfeatures) ? (tid+1)*nfeaturesperthread : nfeatures;
    for(j=pstart; j<pend; j++){
        //fprintf(stderr, "%d-%d\n", cumloci[j], cumloci[j+1]);
        //fprintf(stderr, "%dth peak %d %d\n", j, pstart, pend);
        // eta <- X, beta
#ifdef CBLAS
        fprintf(stderr, "blas\n");
        cblas_dgemv(CblasColMajor, CblasNoTrans, cumloci[j+1]-cumloci[j], P+1, 1.0, X0+cumloci[j], ldx, beta, 1, 0.0, eta+cumloci[j], 1);
#else
        dgemv(X0+cumloci[j], cumloci[j+1]-cumloci[j], P+1, ldx, beta, eta+cumloci[j]);
#endif
        // pjk <- softmax(eta)
        softmax(eta+cumloci[j], pjk+cumloci[j], cumloci[j+1]-cumloci[j]);
        //for(k=cumloci[j]; k<cumloci[j]+10; k++){fprintf(stderr, "%lf ", pjk[k]);}fprintf(stderr, "\n");
        // z   <- pjk, BF, Pi; &
        // Z1  <- z
        tot = 1.0 - Pi[j];
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            z[k] = Pi[j]*pjk[k]*BF[k];
            tot += Pi[j]*pjk[k]*BF[k];
        }
        pmt->lkhd += log(tot);
        //if(isnan(pmt->lkhd)>0 && nanflag==0){fprintf(stderr, "Likelihood became NaN at j=%ld %lf\n", j, Pi[j]); nanflag=1;}
        Z1[j] = 0.0;
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            z[k] /= tot;
            Z1[j] += z[k];
        }
        // w   <- pjk, Z1
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            w[k] = Z1[j] * pjk[k];
        }
        //if(j==5){for(k=cumloci[j]; k<cumloci[j+1]; k++){fprintf(stderr, "%lf %lf ", w[k], Z1[j]);}fprintf(stderr, "\n");}
        // xb  <- X, pjk
#ifdef CBLAS
        cblas_dgemv(CblasColMajor, CblasTrans, cumloci[j+1]-cumloci[j], P, 1.0, X0+cumloci[j], ldx, pjk+cumloci[j], 1, 0.0, xb, 1);
#else
        dgemvT(X0+cumloci[j], cumloci[j+1]-cumloci[j], P, ldx, pjk+cumloci[j], xb);
#endif
        // Xt  <- X, xb, w; &
        // y   <- w, eta, z
        // XjTZj <- X, z
        clearAs0(XjTZj, P);
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            // pseudo data ; P+1th vec of X is offset 
            if(w[k]>0.0){
                y[k] = sqrt(w[k])*(eta[k]-X0[P*ldx+k]) + z[k]/sqrt(w[k]); 
            }else{
                y[k] = 0.0;
            }
            for(l=0; l<P; l++){
                Xt[k+l*(L+P)] = (X0[k+l*ldx]-xb[l]) * sqrt(w[k]);
                XjTZj[l]     += (X0[k+l*ldx]-xb[l]) * z[k];
            }
        }
        // empirical hessian ~= Ic - Io
        for(l=0; l<P; l++){
            for(l2=0; l2<P; l2++){
                Ie[l+l2*P] += XjTZj[l]*XjTZj[l2];
            }
        }
    }
    
    free(XjTZj);
    free(xb);
    
    pthread_exit(NULL);
    return args;
}



int Itr;
void Mstep(double* Xt, double* X0, double* R, double* y, double* beta, double* beta_prior, long L, long P, double lambda){
    long i, j, ldx = P+L;
    for(i=0; i<P; i++){
        y[L+i]=0.0;
        for(j=0; j<P; j++){
            Xt[L+i+ldx*j] = lambda*X0[L+i+ldx*j];
        }
    }
    
    // QR decomposition
    qr(Xt, R, ldx, P, ldx); //fprintf(stderr, "R=\n"); printM2(R, P, P, P);
    // beta1 <- t(Xt) %*% y
#ifdef CBLAS
    cblas_dgemv(CblasColMajor, CblasTrans, ldx, P, 1.0, Xt, ldx, y, 1, 0.0, beta+P+1, 1); fprintf(stderr, "Xt * y=\n"); printV2(beta+P+1, P);
#else
    dgemvT(Xt, ldx, P, ldx, y, beta+P+1); //fprintf(stderr, "Xt * y=\n"); printV2(beta+P+1, P);
#endif
    // beta  <- backsolve(R, beta1)
    BackSolve(R, P, beta+P+1, beta);
    //double th = 15.;
    //for(i=1; i<P; i++){if(beta[i]>th){beta[i]=beta0[i];}else if(beta[i]<(-th)){beta[i]=beta0[i];}}
    for(i=0; i<P; i++){if(isnan(beta[i])>0){fprintf(stderr, "beta[%ld] is NaN\n", i);}}
    if(Itr<1){for(i=0; i<P; i++){beta[i]=beta_prior[i];}}
}





double getQval(double* Z1, double* Pis, double* z, double* pjk, double* BF, long L, long nfeatures, double* X0, double* beta, long Pi, double* U0, double* gamma, long Pj, double* plkhd, double* lambda){
    double pen = 0.0;
    long i, j;
    double tmp;
    long ldx = L+Pi;
    long ldu = nfeatures+Pj;
    // variant level
    for(i=0; i<Pi; i++){
        tmp = 0.0;
        for(j=0; j<Pi; j++){
            tmp += lambda[0]*X0[L+i+j*ldx]*beta[j];
        }
        pen += tmp * tmp;
    }
    // feature level
    for(i=0; i<Pj; i++){
        tmp = 0.0;
        for(j=0; j<Pj; j++){
            tmp += lambda[1]*U0[nfeatures+i+j*ldu]*gamma[j];
        }
        pen += tmp * tmp;
    }
    (*plkhd) += pen;
    return nk_lsum3(Z1, Pis, nfeatures, 1) + nk_lsum2(z, pjk, L, 1) + nk_lsum2(z, BF, L, 1) + pen;
}


long em(double* X0, double* BF, long L, long Pi, double* U0, long nfeatures, long Pj, long* cumloci, double* beta, double* gamma, double* Pis, double* pjk, double* Z1, double* z, long nthreads, double* plkhd){
    long i, j, k, l, l2;
    
    fprintf(stderr, "em: L=%ld Pi=%ld nfeatures=%ld \n", L, Pi, nfeatures);
    
    fprintf(stderr, "Model fitting started...\n\n");
    
    //long Pj = 1;
    //double* U; U = (double*)calloc(nfeatures*Pj, sizeof(double)); for(i=0; i<nfeatures; i++){U[i]=1.;}
    
    
    long tid;
    HIERARCHICAL_MT* pmt; pmt=(HIERARCHICAL_MT*)calloc(nthreads, sizeof(HIERARCHICAL_MT));
    pthread_t* pid;       pid=(pthread_t*)calloc(nthreads, sizeof(pthread_t));
    long nfeaturesperthread  = nfeatures / nthreads + 1;
    
    // parameters
    double sigma = 10.;
    double lambda[2] = {5.0, 1.0};
    //double* Pis;  Pis= (double*)calloc(nfeatures, sizeof(double)); for(i=0; i<nfeatures; i++){Pis[i]=0.075;}
    //double* beta; beta  = (double*)calloc((Pi+1)*2, sizeof(double)); // coef length Pi + 1 + working space for another Pi + 1
    //double* gamma;gamma = (double*)calloc((Pj+1)*2, sizeof(double)); // (Pj + 1) x 2
    double  flip=1.0;
    double* beta_old;  beta_old =  (double*)calloc(Pi+1, sizeof(double)); // coef length Pi + working space for another Pi
    double* gamma_old; gamma_old = (double*)calloc(Pj+1, sizeof(double)); // coef length Pi + working space for another Pi
    beta[Pi] =beta_old[Pi] =beta0[Pi] =1.0;//offset
    gamma[Pj]=gamma_old[Pj]=gamma0[Pj]=1.0;//offset
    cblas_dcopy(Pi, beta0, 1, beta, 1);
    cblas_dcopy(Pi, gamma0, 1, gamma, 1);
    // coproducts
    double* eta;  eta  = (double*)calloc(L, sizeof(double)); // X beta
    //double* Z1;   Z1   = (double*)calloc(nfeatures, sizeof(double)); // 1-Z0j; j=1,...,nfeatures
    //double* z;    z    = (double*)calloc(L, sizeof(double)); // porsterior probs
    double* y;    y    = (double*)calloc(L, sizeof(double)); // pseudo data : W^1/2 %*% eta + W^1/2 %*% z
    double* w;    w    = (double*)calloc(L, sizeof(double)); // weights for IRLS : expand(Z1) * pjk
    //double* pjk;  pjk  = (double*)calloc(L, sizeof(double));  // softmax prior
    double* Xt;   Xt   = (double*)calloc((L+Pi)*Pi, sizeof(double)); // normalized X & replaced by Q
    double* R;    R    = (double*)calloc(Pi*Pi, sizeof(double)); // upper tri
    double* Ie;   Ie   = (double*)calloc(Pi*Pi*nthreads, sizeof(double));
    fprintf(stderr, "Memory allocated...\n\n");
    
    // for peak
    double* Ut;   Ut   = (double*)calloc((nfeatures+Pj)*(Pj+1), sizeof(double)); // model matrix (W^1/2 U, Phi^-1/2)
    double* y2;   y2   = (double*)calloc( nfeatures+Pj   , sizeof(double)); // pseudo data
    double* R2;   R2   = (double*)calloc(Pj*Pj, sizeof(double)); // upper tri
    double* IeU;  IeU  = (double*)calloc(Pj*Pj, sizeof(double));
    
    
    
    
    
    
    long itr;
    
    int again=0;
    double lkhd, lkhd1=-1.0e10;
    double qval, qval1=-1.0e10;
    int conv=0;
    for(itr=0; itr<1000; itr++){
        Itr=itr;
        //############
        //#  Esteps  #
        //############
        for(tid=0; tid<nthreads; tid++){
            pmt[tid].tid    = tid;
            pmt[tid].nfeaturesperthread = nfeaturesperthread;
            pmt[tid].nfeatures = nfeatures;
            pmt[tid].cumloci = cumloci;
            pmt[tid].Pi = Pi;
            pmt[tid].X0 = X0;
            pmt[tid].BF = BF;
            pmt[tid].beta = beta;
            pmt[tid].eta  = eta;
            pmt[tid].Pis  = Pis;
            pmt[tid].pjk  = pjk;
            pmt[tid].z    = z;
            pmt[tid].Z1   = Z1;
            pmt[tid].w    = w;
            pmt[tid].y    = y;
            pmt[tid].Xt   = Xt;
            pmt[tid].Ie   = Ie+Pi*Pi*tid;
            pmt[tid].lkhd = 0.0;
            //Estep(pmt);
            
            long pthflag;
            if( (pthflag = pthread_create(pid+tid, NULL, (void*)Estep, (void*)(pmt+tid))) != 0){
                fprintf(stderr, "Thread not created...aborted.\n");
                return -9999;
            }
        }
        for(tid=0; tid<nthreads; tid++){
            long pthflag;
            if( (pthflag = pthread_join(pid[tid], NULL)) !=0 ){fprintf(stderr, "Thread not joined...aborted.\n"); return -9999;};
        }
        
        
        lkhd = 0.0;
        for(tid=0; tid<nthreads; tid++){lkhd += pmt[tid].lkhd; if(isnan(pmt[tid].lkhd)>0){fprintf(stderr, "Partial likelihood for thread %d is NaN.\n", tid);} }
        qval = getQval(Z1, Pis, z, pjk, BF, L, nfeatures, X0, beta, Pi, U0, gamma, Pj, &lkhd, lambda);
        //fprintf(stderr, "qval=%lf -> %lf\n", qval1, qval);
        fprintf(stderr, "lkhd=%lf -> %lf\n", lkhd1, lkhd);
        fprintf(stderr, "Current : "); printV2(beta, Pi+1); printV2(gamma, Pj+1);
        if(isnan(lkhd)>0 || (lkhd1-lkhd>100.0 && again<20)){
            //fprintf(stderr, "  devol!\n");
            fprintf(stderr, "  Before: "); printV2(beta, Pi);
            for(i=0; i<Pi; i++){ beta[i]  = beta_old[i]  + flip*(beta[i]  - beta_old[i]) /(flip>0.0 ? 10.0 : 1.0); }
            for(i=0; i<Pj; i++){ gamma[i] = gamma_old[i] + flip*(gamma[i] - gamma_old[i])/(flip>0.0 ? 10.0 : 1.0); }
            dgemv(U0, nfeatures, Pj+1, nfeatures+Pj, gamma, Pis); for(i=0; i<nfeatures; i++){Pis[i]=1./(1.+exp(-Pis[i]));}
            fprintf(stderr, "  After : "); printV2(beta, Pi);
            flip *= (-1.0);
            again++;
        }else{
            //###############
            //# Mstep for X #
            //###############
            cblas_dcopy(Pi, beta,  1, beta_old,  1);
            cblas_dcopy(Pj, gamma, 1, gamma_old, 1);
            fprintf(stderr, "Mstep\n");
            if(Pi>0) Mstep(Xt, X0, R, y, beta, beta0, L, Pi, lambda[0]);
            fprintf(stderr, "beta="); printV2(beta, Pi);
            for(tid=1; tid<nthreads; tid++){for(l=0; l<Pi; l++){for(l2=0; l2<Pi; l2++){Ie[l+l2*Pi]+=Ie[l+l2*Pi+tid*Pi*Pi];}}}
            
            //###############
            //# Mstep for U #
            //###############
            for(i=0; i<nfeatures; i++){
                y2[i] = (Z1[i]-Pis[i])/sqrt(Pis[i]*(1.-Pis[i]));
                for(j=0; j<Pj; j++){
                    Ut[i+j*(nfeatures+Pj)] = U0[i+j*(nfeatures+Pj)] * sqrt(Pis[i]*(1.-Pis[i]));
                    y2[i] += Ut[i+j*(nfeatures+Pj)]*gamma[j];
                }
            }
            Mstep(Ut, U0, R2, y2, gamma, gamma0, nfeatures, Pj, lambda[1]);
            fprintf(stderr, "gamma=");
            printV2(gamma, Pj);
            dgemv(U0, nfeatures, Pj+1, nfeatures+Pj, gamma, Pis);// Pi <- U %*% gamma // not yet Pi
            for(i=0; i<nfeatures; i++){Pis[i]=1./(1.+exp(-Pis[i]));} // Pis <- 1/(1+exp(Pis)))
            
            // terminate
            if(fabs(lkhd-lkhd1)<1.0e-4 && again==0 && itr>10){conv=1; break;}
            flip = 1.0;
            lkhd1=lkhd;
            again = 0;
        }
    }
    
    //##############
    //# Output
    //##############
    
    (*plkhd) = lkhd;

    // hessian for beta
    for(j=0; j<Pi*Pi; j++){beta[Pi+j]=Ie[j];}
    // hessian for gamma
    clearAs0(IeU, Pj*Pj);
    for(i=0; i<nfeatures; i++){
        double zmp = (Z1[i]-Pis[i]);
        for(l=0; l<Pj; l++){
            for(l2=0; l2<Pj; l2++){
                IeU[l+l2*Pj] += U0[i+l*(nfeatures+Pj)]*U0[i+l2*(nfeatures+Pj)]*zmp*zmp;
            }
        }
    }
    for(j=0; j<Pj*Pj; j++){gamma[Pj+j]=IeU[j];}
    
    return (conv==1 ? itr : -itr);

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

long main(long argc, char** argv){
    long i, j, k, l;
   
    if(argc==1){usage_hm(); return 1;} 
    long verbose = 0; for(i=0; i<argc; i++){if(strcmp(argv[i], "-v")==0){ verbose=1; }}
    
    gzFile fi = NULL; // i = variant level
    gzFile fj = NULL; // j = feature level
    char* typi = NULL;
    char* typj = NULL;
    long* nxpi = NULL;
    long* nxpj = NULL;
    long* cumcoli;
    long* cumcolj;
    long* cumloci;
    long nvari, nvarj;
    long nfeatures=0;
    long pi, pj, Pi, Pj;// numbers of col and expanded col (p & P)
    long nrow=2;
    
    double** Xki;
    double** Xkj;
    double** Bsi;
    double** Bsj;
    
    double* X;
    double* U;
    double* bf;
    
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-f")==0){ nfeatures = (long)atoi(argv[i+1]); }}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-r")==0){ nrow = (long)atoi(argv[i+1]); }}
    
    
    // input for variant level
    for(i=0; i<argc-1; i++){
        if(strcmp(argv[i], "-i")==0){
            fi = gzopen(argv[i+1], "rb6f");
            for(j=0; j<argc-1; j++){
                if(strcmp(argv[j], "-c")==0){
                    pi = parseCol(argv[j+1], &typi, &nxpi);
                    Pi = countP(&typi, &nxpi, pi);
                    Xki = (double**)calloc(pi, sizeof(double*));
                    Bsi = (double**)calloc(pi, sizeof(double*));
                    cumcoli = (long*)calloc(Pi+1, sizeof(long));cumcoli[0]=0;
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
                        }
                    }
                    break;
                }
            }
            nvari = l;
            if(typi==NULL){fprintf(stderr, "No column information for %s\n", argv[i+1]); return 1;}
            
            X = (double*)calloc((nrow+Pi)*(Pi+1), sizeof(double));
            bf= (double*)calloc(nrow,        sizeof(double));
            cumloci = (long*)calloc(nfeatures+1, sizeof(double));
            readTable(fi, X, bf, typi, nxpi, Xki, Bsi, nrow, pi, cumloci, cumcoli, Pi, 100.0);
            if(verbose>0){ 
                //fprintf(stderr, "cumloci: "); printVL(cumloci, nfeatures+1);
                //fprintf(stderr, "cumcoli: "); printVL(cumcoli, nvari+1);
                //printM(X, nrow+Pi, Pi+1);
                printM2(X+nrow, Pi, Pi+nrow, Pi);
                //fprintf(stderr, "BF: "); printV(bf, 10);
            }
            break;
        }
    }
    if(fi==NULL){fprintf(stderr, "No input file!\n"); return 1;}
    
    // input for feature level
    for(i=0; i<argc-1; i++){
        if(strcmp(argv[i], "-j")==0){
            fj = gzopen(argv[i+1], "rb6f");
            for(j=0; j<argc-1; j++){
                if(strcmp(argv[j], "-d")==0){
                    pj = parseCol(argv[j+1], &typj, &nxpj);
                    Pj = countP(&typj, &nxpj, pj)+1;// +1 for intersept
                    Xkj = (double**)calloc(pj, sizeof(double*));
                    Bsj = (double**)calloc(pj, sizeof(double*));
                    cumcolj = (long*)calloc(Pj+1, sizeof(long)); cumcolj[0]=1;
                    l=0;
                    for(k=0; k<pj; k++){
                        if(verbose>0)fprintf(stderr, "%ld %c %ld\n", Pj, typj[k], nxpj[k]);
                        if(typj[k]=='N' || typj[k]=='C'){
                            cumcolj[l+1] = cumcolj[l] + (typj[k]=='N' ? nxpj[k]+1 : nxpj[k]-1);
                            if(typj[k]=='N' && nxpj[k]>0){
                                getXk(&(Xkj[l]),nxpj[k]); 
                                getBS(&(Bsj[l]), Xkj[l], nxpj[k]);
                                if(verbose>0){fprintf(stderr, "knots: "); printV2(Xkj[l], nxpj[k]);}
                            }
                            if(verbose>0) fprintf(stderr, "cumcolj=%ld\n", cumcolj[l+1]);
                            l++;
                        }
                    }
                    break;
                }
            }
            nvarj = l;
            if(typj==NULL){fprintf(stderr, "No column information for %s\n", argv[i+1]); return 1;}
            
            U = (double*)calloc((nfeatures+Pj)*(Pj+1), sizeof(double)); for(j=0; j<nfeatures; j++){U[j]=1.0;}
            readTable(fj, U, NULL, typj, nxpj, Xkj, Bsj, nfeatures, pj, NULL, cumcolj, Pj, 100.0);
            if(verbose>0){
                //printM(U, nfeatures+Pj, Pj+1);
                printM2(U+nfeatures, Pj, Pj+nfeatures, Pj);
                //fprintf(stderr, "cumcolj: "); printVL(cumcolj, nvarj+1);
            }
            break;
        }
    }
    if(fj==NULL){
        fprintf(stderr, "No covariate for feature level prior probability.\n"); 
        Pj = 1;
        U = (double*)calloc((nfeatures+Pj)*(Pj+1), sizeof(double)); for(j=0; j<nfeatures; j++){U[j]=1.0;}
    }
    
    double* beta; beta  = (double*)calloc((Pi+1)*2+Pi*Pi, sizeof(double)); // coef length Pi + 1 + working space for another Pi + 1
    double* gamma;gamma = (double*)calloc((Pj+1)*2+Pj*Pj, sizeof(double)); // (Pj + 1) x 2
    double* Z1;   Z1   = (double*)calloc(nfeatures, sizeof(double)); // 1-Z0j; j=1,...,nfeatures
    double* z;    z    = (double*)calloc(nrow, sizeof(double)); // porsterior probs
    double* pjk;  pjk  = (double*)calloc(nrow, sizeof(double));  // softmax prior
    
    double* Pis;  Pis= (double*)calloc(nfeatures, sizeof(double)); for(i=0; i<nfeatures; i++){Pis[i]=0.075;}
    
    // Init parameter setting
    beta0 =  (double*)calloc((Pi+1)*2, sizeof(double)); // coef length Pi + 1 + working space for another Pi + 1
    gamma0 = (double*)calloc((Pj+1)*2, sizeof(double)); // (Pj + 1) x 2
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--init-variant-level")==0){
        FILE* fbeta0;
        fbeta0 = fopen(argv[i+1], "rb");
        fread(beta0, sizeof(double), Pi, fbeta0);
        fclose(fbeta0);
        fprintf(stderr, "Init Beta: "); printV2(beta0, Pi);
    }}
    //beta0[0] = -0.768617;
    //beta0[1] =  3.072989;
    //if(Pi>2)beta0[2] = 2.631762;
    
    gamma0[0] = -2.6;
    /*gamma0[0] = -8.000779;
    gamma0[1] = 9.790940;
    gamma0[2] = 145.083555;
    gamma0[3] = 95.627158;
    gamma0[4] = 335.465419;
    gamma0[5] = 206.714531;*/
    
    
    
    
    
    long nthreads = 1;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-t")==0){ nthreads = (long)atoi(argv[i+1]); }}
    
    double lkhd;
    long itr;
    itr = em(X, bf, nrow, Pi, U, nfeatures, Pj, cumloci, beta, gamma, Pis, pjk, Z1, z, nthreads, &lkhd);
    
    
    
    
    
    
    
    
    
    
    
    //##########
    //# Output #
    //##########
    
    
    char* prefix;
    char* ofname;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-o")==0){ prefix = argv[i+1]; }}
    ofname = (char*)calloc(strlen(prefix)+100, sizeof(char));
    
    // directry exists?
    struct stat st = {0}; if (stat(prefix, &st) == -1) { mkdir(prefix, 0700); }
    
    // log file
    sprintf(ofname, "%s/log.txt", prefix);
    FILE* flog;  flog  = fopen(ofname, "wb");
    fprintf(flog, "Likelihood   : %lf\n", lkhd);
    fprintf(flog, "N iterations : %ld\n", itr);
    fprintf(flog, "N threads    : %ld\n", nthreads);
    fprintf(flog, "Variant-level: "); for(i=0; i<Pi-1; i++){fprintf(flog, "%lf, ", beta[i]);} fprintf(flog, "%lf\n", beta[i]);
    fprintf(flog, "Feature-level: "); for(i=0; i<Pj-1; i++){fprintf(flog, "%lf, ", gamma[i]);} fprintf(flog, "%lf\n", gamma[i]);

    fclose(flog);
    
    sprintf(ofname, "%s/variant_level.bin", prefix);
    FILE* hessf;  hessf  = fopen(ofname, "wb");
    sprintf(ofname, "%s/feature_level.bin", prefix);
    FILE* hessf2; hessf2 = fopen(ofname, "wb");
    sprintf(ofname, "%s/Pi1.bin", prefix);
    FILE* pibin;  pibin  = fopen(ofname,   "wb");
    fwrite(beta,  Pi*(Pi+1), sizeof(double), hessf);
    fwrite(Pis,   nfeatures, sizeof(double), pibin);
    fwrite(gamma, Pj*(Pj+1), sizeof(double), hessf2);
    
    
    
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-p")==0){
        sprintf(ofname, "%s/pp.gz", prefix);
        gzFile postf; postf = gzopen(ofname, "wb6f");
        sprintf(ofname, "%s/Z1.gz", prefix);
        gzFile zf;    zf    = gzopen(ofname, "wb6f");
        for(j=0; j<nfeatures; j++){
            gzprintf(zf, "%lf\n", Z1[j]);
            for(k=cumloci[j]; k<cumloci[j+1]; k++){
                gzprintf(postf, "%ld\t%lf\t%lf\t%lf\n", j+1, pjk[k], bf[k], z[k]);
            }
        }
        gzclose(postf);
        gzclose(zf);
    }}
    // number of loci that acounts for postth% of posterior prob
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--credible-set")==0){
        sprintf(ofname, "%s/credible_set.gz", prefix);
        gzFile nplocif; nplocif = gzopen(ofname, "wb6f");
        double postth[]={0.9, 0.95, 0.99};
        int    nlocij[]={0,0,0};
        double plocij=0.0;
        for(j=0; j<nfeatures; j++){
            qsort(z+cumloci[j], cumloci[j+1]-cumloci[j], sizeof(double), cmpdr);
            plocij = 0.0;
            nlocij[0] = nlocij[1] = nlocij[2] = 0;
            for(k=cumloci[j]; k<cumloci[j+1]; k++){
                //if(j==1){fprintf(stderr, "%lf ", z[k]/Z1[j]);}
                if(plocij<postth[0]){nlocij[0]++;}
                if(plocij<postth[1]){nlocij[1]++;}
                if(plocij<postth[2]){nlocij[2]++;}
                plocij += z[k]/Z1[j];
            }
            gzprintf(nplocif, "%ld\t%d\t%d\t%d\n", j+1, nlocij[0], nlocij[1], nlocij[2]);
        }
        gzclose(nplocif);
    }}
    
}









