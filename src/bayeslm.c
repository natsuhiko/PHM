#include "bayeslm.h"
#include "util_bayeslm.h"
#include "usage.h"


int startWith(const char *pre, const char *str)
{
    return strncmp(pre, str, strlen(pre));
}

void clearAs(double* x, int n, double val){
    int i;
    for(i=0; i<n; i++){
        x[i] = val;
    }
}



double nk_var(double* x, double* y, int N){
    double mx, my;
    int i;
    mx=my=0.0;
    double res=0.0;
    for(i=0; i<N; i++){
        mx += x[i]/(double)N;
        my += y[i]/(double)N;
    }
    for(i=0; i<N; i++){
        res += (x[i]-mx)*(y[i]-my);
    }
    return res/(double)N;
};


double nk_cor(double* x, double* y, int n){
    return nk_var(x, y ,n) / sqrt(nk_var(x, x, n)) / sqrt(nk_var(y, y, n));
}




double nk_lm2(double* X, int* id, double* y, int n, double* work){
    clearAs0(work, 9);
    double* xtx = work;
    double* xy  = work + 3;
    double* ms  = work + 6;
    int i=0, j;
    double dn = (double)n;
    
    // mean
    for(i=0; i<n; i++){
        ms[0] += X[id[0]*n+i] / dn;
        ms[1] += X[id[1]*n+i] / dn;
        ms[2] += y[i] / dn;
    }
    // var
    for(i=0; i<n; i++){
        xtx[0] += (X[id[0]*n+i]-ms[0])*(X[id[0]*n+i]-ms[0]);
        xtx[1] += (X[id[0]*n+i]-ms[0])*(X[id[1]*n+i]-ms[1]);
        xtx[2] += (X[id[1]*n+i]-ms[1])*(X[id[1]*n+i]-ms[1]);
        
        xy[0]  += (X[id[0]*n+i]-ms[0])*(y[i]-ms[2]);
        xy[1]  += (X[id[1]*n+i]-ms[1])*(y[i]-ms[2]);
        xy[2]  += (y[i]-ms[2])        *(y[i]-ms[2]);
    }
    // inverse
    double det = xtx[0]*xtx[2]-xtx[1]*xtx[1];
    if(det<=0.0){return 0.0;}
    double tmp =  xtx[2]/det;
    xtx[2]     =  xtx[0]/det;
    xtx[0]     =  tmp;
    xtx[1]     = -xtx[1]/det;
    
    double res = 0.0;
    for(i=0; i<2; i++){
        for(j=0; j<2; j++){
            res += xy[i]*xy[j]*xtx[i+j];
        }
    }
    return sqrt(res/xy[2]);
}

int parseFormat(char* str, int* formatID){
    int i;
    int nfield=1;
    for(i=0; i<strlen(str); i++){if(str[i]==':'){nfield++;}}
    char format[3];
    
    int ns=0;
    for(i=0; i<nfield; i++){
        if(i<nfield-1){
            sscanf(str+ns, "%[^:]:", format);
        }else{
            sscanf(str+ns, "%s", format);
        }
        ns += strlen(format)+1;
        //fprintf(stderr, "%d %s ", ns, format);
        if(strcmp(format,"GT")==0){
            formatID[i]=FORMAT_GT;
        }else if(strcmp(format,"GL")==0){
            formatID[i]=FORMAT_GL;
        }else if(strcmp(format,"AP")==0){
            formatID[i]=FORMAT_AP;
        }else if(strcmp(format,"GP")==0){
            formatID[i]=FORMAT_GP;
        }else if(strcmp(format,"PP")==0){
            formatID[i]=FORMAT_GP;
        }else if(strcmp(format,"AS")==0){
            formatID[i]=FORMAT_AS;
        }else if(strcmp(format,"RD")==0){
            formatID[i]=FORMAT_RD;
        }else if(strcmp(format,"BF")==0){
            formatID[i]=FORMAT_BF;
        }else if(strcmp(format,"DS")==0){
            formatID[i]=FORMAT_DS;
        }else{
            formatID[i]=FORMAT_OTHER;
        }
    }
    return nfield;
}

int doseFormatExist(int* formatID, int nfield, int formatID1){
    int i;
    for(i=0; i<nfield; i++){
        if(formatID[i]==formatID1){return 1;}
    }
    return 0;
}

int parseInfo(char* str, char* infostr, VCF_info* vinfo){
    int i;
    int nfield=1;
    for(i=0; i<strlen(str); i++){if(str[i]==';'){nfield++;}}
    int ns=0;
    vinfo->VT=-100;
    vinfo->RSQ = -1.0;
    for(i=0; i<nfield; i++){
        sscanf(str+ns, "%[^;];", infostr);
        if(strcmp(infostr,"VT=SNP")==0){
            vinfo->VT=VT_SNP;
        }else if(strcmp(infostr,"VT=INDEL")==0){
            vinfo->VT=VT_INDEL;
        }else if(strcmp(infostr,"VT=SV")==0){
            vinfo->VT=VT_SV;
        }else if(startWith("RSQ=",infostr)==0){
            sscanf(infostr, "RSQ=%lf", &(vinfo->RSQ));
        }else if(startWith("AF=",infostr)==0){
            sscanf(infostr, "AF=%lf", &(vinfo->AF));
        }else if(startWith("IMP2=",infostr)==0){
            sscanf(infostr, "IMP2=%lf,%lf,%lf", &(vinfo->AF), &(vinfo->RSQ), &(vinfo->CER));
        }
        ns += strlen(infostr)+1;
    }
    return 0;
}



int isSnp(int* allen, int nal){
    int i;
    for(i=0; i<nal; i++){
        if(allen[i]>1){return 0;}
    }
    return 1;
}


/*
double nk_dsum(double* x, int n, int ldx){
    int i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}
*/


int nk_isum(int* x, int n, int ldx){
    int i;
    int res=0;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

void gt2ap(int* gt, int n, double* ap){
    int i;
    double* p; 
    double* q;
    clearAs0(ap, 2*n);
    p = ap;
    q = ap+n;
    p[gt[0]] = 1.0;
    q[gt[1]] = 1.0;
}

void rsq2ap(int* gt, double rsq, double* afs, int n, int samplesize, double* ap){// all samples with afs
    int i;
    double tafs = 0.0;
    for(i=0; i<n; i++){
        tafs += afs[i]*(1.0-afs[i]);
    }
    double ep = (1.0-rsq) * tafs / (2.0*((double)(n-1)));
    double* p; 
    double* q;
    clearAs(ap, 2*n*samplesize, ep);
    for(i=0; i<samplesize; i++){
        p = ap+i*n*2;
        q = ap+i*n*2+n;
        p[gt[i*2+0]] = 1.0-ep*(double)(n-1);
        q[gt[i*2+1]] = 1.0-ep*(double)(n-1);
    }
}

void ds2ap(int* gt, double* ds, int n, double* ap){
    if(gt[0]==gt[1]){
        int k;
        for(k=0; k<n; k++){
            ap[k]=ap[k+n]=ds[k]/2.0;
        }
    }else{
        clearAs0(ap, n*2);
        double* p; p = ap;
        double* q; q = ap+n;
        p[gt[0]] = ds[gt[0]];
        q[gt[1]] = ds[gt[1]];
        if(p[gt[0]]>1.0){p[gt[0]]=1.0; q[gt[0]]=ds[gt[0]]-1.0;}
        if(q[gt[1]]>1.0){q[gt[1]]=1.0; p[gt[1]]=ds[gt[1]]-1.0;}
        
        double tp = 1.0 - p[gt[0]] - p[gt[1]];
        double tq = 1.0 - q[gt[0]] - q[gt[1]];
        
        int k;
        for(k=0; k<n; k++){
            if(k!=gt[0] && k!=gt[1]){
                p[k] = tp/(tp+tq)*ds[k];
                q[k] = tq/(tp+tq)*ds[k];
            }
        }
    }
}

void gl2ap(int* gt, double* gl, int n, double* ap, double* d){
    clearAs(ap, n*2, 0.1/((double)(n-1)));
    double* p; p = ap;
    double* q; q = ap+n;
    p[gt[0]] = 0.9;
    q[gt[1]] = 0.9;
    int itr, j, k, l;
    double denom, lkhd, lkhd0=-1.0e10;
    for(itr=0; itr<100; itr++){
        clearAs0(d, n*n);
        l=0;
        lkhd=0.0;
        for(k=0; k<n; k++){
            for(j=0; j<=k; j++){
                denom = p[j]*q[k] + p[k]*q[j];
                if(denom>1.0e-20){
                    lkhd     += gl[l]*log(denom);
                    d[j*n+k] += gl[l]*p[j]*q[k] / denom;
                    d[k*n+j] += gl[l]*p[k]*q[j] / denom;
                }
                l++;
            }
        }
        for(j=0; j<n; j++){
            p[j] = nk_dsum(d+j*n, n, 1);
            q[j] = nk_dsum(d+j,   n, n);
        }
        if(fabs(lkhd-lkhd0)<1.0e-8){ break; }else{ lkhd0=lkhd; }
    }
}

int choose(int n, int k){
    int i;
    int res = 1;
    for(i=n; i>=n-k+1; i--){
        res *= i;
    }
    for(i=k; i>=1; i--){
        res /= i;
    }
    return res;
}
int achoose(int n){
    int i;
    int res=0;
    for(i=1; i<=n-1; i++){
        res += choose(n, i);
    }
    return res/2;
}

int getCombAlk(int n, int k, int* v, int id, int idmax, char* a0, char* a1, char** ba0, char** ba1){
    if(k==n){
        if(id==-1 || id>=idmax){return id+1;}
        int i, i0, i1, allen=0;
        memcpy(ba0[id], a0, strlen(a0));
        i0=strlen(a0);
        i1=0;
        k=1; // kth allele for a1
        for(i=0; i<strlen(a1)+1; i++){
            if(a1[i]==',' || a1[i]=='\0'){
                if(v[k]==1){
                    if(i1==0){
                        memcpy(ba1[id]+i1, a1+(i-allen), allen);
                        i1 += allen;
                    }else{
                        ba1[id][i1++]=',';
                        memcpy(ba1[id]+i1, a1+(i-allen), allen);
                        i1 += allen;
                    }
                }else{
                    ba0[id][i0++]=',';
                    memcpy(ba0[id]+i0, a1+(i-allen), allen);
                    i0 += allen;
                }
                k++;
                allen=0;
            }else{
                allen++;
            }
        }
        ba0[id][i0]='\0';
        ba1[id][i1]='\0';
        
        return id+1;
    }
    v[k]=0;
    id = getCombAlk(n, k+1, v, id, idmax, a0, a1, ba0, ba1);
    v[k]=1;
    id = getCombAlk(n, k+1, v, id, idmax, a0, a1, ba0, ba1);
}


int getCombk(int n, int k, int* v, int id, int idmax, double* ds, double* bds, int samplesize){
    if(k==n){
        if(id==-1 || id>=idmax){return id+1;}
        int i;
        //printf("%d ", id); for(i=0; i<n; i++){printf("%d ", v[i]);} 
        for(i=0; i<n; i++){
            if(v[i]>0){
                bds[id*samplesize]+=ds[i];
            }
        }
        //printf("%lf\n", bds[id]);
        
        return id+1;
    }
    v[k]=0;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize);
    v[k]=1;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize);
}


int countFields(char* alStr, char sep){
    int i;
    int n=0;
    for(i=0; i<strlen(alStr); i++){
        if(alStr[i]==sep){
            n++;
        }
    }
    return n+1;
}

int isPhased(char* x, int l){
    int i;
    for(i=0; i<l; i++){
        //fprintf(stderr, "%c\n", x[i]);
        if(x[i]=='|'){
            return 1;
        }else if(x[i]=='/'){
            return 0;
        }
    }
}

int parseGT(char* cell, int* gt){
    if(cell[0]=='.' && cell[1]=='/' && cell[2]=='.'){gt[0]=gt[1]=-9999999; return 4;}
    if(cell[0]=='.'){gt[0]=gt[1]=-9999999; return 2;}
    int nchar1;
    if(isPhased(cell, strlen(cell))>0){
        sscanf(cell, "%d|%d%n", gt, gt+1, &nchar1);
    }else{
        sscanf(cell, "%d/%d%n", gt, gt+1, &nchar1);
    }
    return nchar1 + 1;
}

int parseDS(char* cell, int n, double* ds){
    if(cell[0]=='.'){ds[0]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=1; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", ds+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", ds+i, &nchar1);
    
    ds[0] = 1.0 - nk_dsum(ds+1, n-1, 1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int skipFormat(char* cell){
    int l=strlen(cell);
    int i;
    for(i=0; i<l; i++){
        if(cell[i]=='\t' || cell[i]==':' || cell[i]=='\n' || cell[i]=='\0'){return i+1;}
    }
}

int parseGP(char* cell, int n, double* gl){
    if(cell[0]=='.'){gl[0]=gl[1]=gl[2]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(cell+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1 + 1;
    return nchar;
}

int parseGL(char* cell, int n, double* gl){
    if(cell[0]=='.'){gl[0]=gl[1]=gl[2]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(cell+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1 + 1;
    for(i=0; i<n*(n-1)/2+n; i++){
        gl[i] = pow(10.0, gl[i]);
    }
    return nchar;
}

int parseAP(char* cell, int n, double* ap){
    if(cell[0]=='.'){ap[0]=ap[1]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=1; i<n; i++){
        sscanf(cell+nchar, "%lf,%n", ap+i, &nchar1);
        nchar += nchar1;
    }
    for(i=n+1; i<2*n-1; i++){
        sscanf(cell+nchar, "%lf,%n", ap+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", ap+i, &nchar1);
    
    ap[0] = 1.0 - nk_dsum(ap+1,   n-1, 1);
    ap[n] = 1.0 - nk_dsum(ap+n+1, n-1, 1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int parseAS(char* cell, int n, double* as){
    if(cell[0]=='.'){as[0]=as[1]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=0; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", as+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", as+i, &nchar1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int parseCell(char* cell, int nal, int* gt, double* gl, double* ap, double* asc, double* ds, int* formatID, int nfields){
    int i, offs=0, k=0;
    for(k=0; k<nfields; k++){
        if(formatID[k]==FORMAT_GT){
            offs += parseGT(cell+offs, gt);
        }else if(formatID[k]==FORMAT_GL){
            offs += parseGL(cell+offs, nal, gl);
        }else if(formatID[k]==FORMAT_GP){
            offs += parseGP(cell+offs, nal, gl);
        }else if(formatID[k]==FORMAT_AP){
            offs += parseAP(cell+offs, nal, ap);
        }else if(formatID[k]==FORMAT_AS){
            offs += parseAS(cell+offs, nal, asc);
        }else if(formatID[k]==FORMAT_DS){
            offs += parseDS(cell+offs, nal, ds);
        }else{
            offs += skipFormat(cell+offs);
        }
    }
    return offs;
}




int parseBody3(char* body, int n, int* gt, double* ds, double* gl, double* ap, double* d){
    int i;
    int nchar=0, nchar1;
    
    // GT
    sscanf(body+nchar, "%d|%d:%n", gt, gt+1, &nchar1);
    nchar += nchar1;
    
    // DS
    for(i=1; i<n-1; i++){
        sscanf(body+nchar, "%lf,%n", ds+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%lf:%n", ds+i, &nchar1);
    nchar += nchar1;
    
    // GL
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(body+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1;
    
    return nchar+1;
}

int parseBody1(char* body, int n, int* gt, double* ds, double* gl, double* ap, double* d){
    int nchar1;
    sscanf(body, "%d|%d%n", gt, gt+1, &nchar1);
    return nchar1+1;
}

int parseBody(char* body, int n, int* gt){
    int i;
    int nchar=0, nchar1;
    for(i=0; i<n-1; i++){
        sscanf(body+nchar, "%d|%d\t%n", gt+2*i, gt+2*i+1, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%d|%d", gt+2*i, gt+2*i+1);
    return i+1;
}

void gt2dsgl(int* gt, int nal, double* ds, double* gl){
    int j, k, l;
    clearAs0(ds, nal);
    clearAs0(gl, choose(nal,2)+nal);
    ds[gt[0]]++;
    ds[gt[1]]++;
    l=0;
    for(k=0; k<nal; k++){
        for(j=0; j<=k; j++){
            if((j==gt[0] && k==gt[1]) || (k==gt[0] && j==gt[1])){gl[l]++; break;}
            l++;
        }
    }
}

int parseCigarPoint(int start, char* c1, char* seq, int* at, int K, char*** als, int** asc, int* nals, int** allen){
    //char* op; // operation
    //op = (char*)calloc(100, sizeof(char));
    char op[10];
    int len;// length in cigar
    int start0=start;
    int nseq=0;
    char* pc1;
    int nc1=0;
	int ali, ai, k;
    //pc1 = c1;
    int end;
	int flag=0;
	int flagin=0;
	int nchar=0;
    int allelematch;
    char prevop;
    while((sscanf(c1+nc1, "%d%[MIDNSH=PX]%n", &len, op, &nchar))>0){
        nc1 += nchar;
		flagin=0;
        if(op[0]=='M' || op[0]=='=' || op[0]=='X'){
            end = start+len-1;
			int nseq0=nseq;
			for(k=0; k<K; k++){
				nseq=nseq0;
                if(start <= at[k] && at[k] <= end){
					flagin++;
                    nseq += (at[k] - start);
                    for(ai=0; ai<nals[k]; ai++){
#ifdef INDEL
                        if(end-at[k]+1 >= allen[k][ai]){
                            allelematch=1;
                            for(ali=0; ali<allen[k][ai]; ali++){
                                if(seq[nseq+ali]!=als[k][ai][ali]){
                                    allelematch=0;
                                    break;
                                }
                            }
                        }else{// allele exceeds match region
                            allelematch=0;
                        }
                        asc[k][ai] += allelematch;
#else
                        if(seq[nseq]==als[k][ai][0] && allen[k][ai]==1){
                            asc[k][ai]++;
                        }
#endif
                    }
                }
			}
            start += len;
            nseq = nseq0+len;
        }else if(op[0]=='D'){
#ifdef INDEL
            end = start+len-1;
            for(k=0; k<K; k++){
                if(start-1 == at[k] && end < at[k]+allen[k][0]){
                    fprintf(stderr, "%s\n", seq+nseq-1);
					flagin++;
                    for(ai=1; ai<nals[k]; ai++){// only alt allele(s)
                        if(allen[k][0]-allen[k][ai]==len){
                            //fprintf(stderr, "%d %d\n", allen[k][0], allen[k][ai]);
                            asc[k][ai]++;
                        }
                    }
                }
			}
#endif
            start += len;
        }else if(op[0]=='N' || op[0]=='P'){
            start += len;
        }else{// S, H and I
            nseq += len;
        }
		if(flagin>0){flag++;}
        prevop = op[0];
    }
    
    return flag;
}




void createDataBase(const char* fname, const char* reg){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nrs;
    char* rs1; rs1 = (char*)calloc(1000, sizeof(char));
    int rspos;
    char* rschr; rschr = (char*)calloc(1000, sizeof(char));
    int geta;
    char* num; num = (char*)calloc(100, sizeof(char));
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        nrs = countFields(rs, ';');
        if(nrs==1){
            if(strlen(rs)>3){
                sscanf(rs, "%[^0-9]%d", rschr, &rspos);
                geta=0;
                if(rspos>500000000){geta = sprintf(num, "%d", rspos/100000000); rspos = rspos%100000000;}
                rs[4+geta] = '\0';
                fprintf(stdout, "%s\t%d\t%s\t%d\t%s\t%s\n", rs, rspos, chr, pos, a0, a1);
            }
        }else{
            int rslen=0;
            for(i=0; i<strlen(rs)+1; i++){
                if(rs[i]==';' || rs[i]=='\0'){
                    if(rslen>3){
                        memcpy(rs1, rs+(i-rslen), rslen);
                        rs1[rslen]='\0';
                        sscanf(rs1, "%[^0-9]%d", rschr, &rspos);
                        geta=0;
                        if(rspos>500000000){geta = sprintf(num, "%d", rspos/100000000); rspos = rspos%100000000;}
                        rs1[4+geta] = '\0';
                        fprintf(stdout, "%s\t%d\t%s\t%d\t%s\t%s\n", rs1, rspos, chr, pos, a0, a1);
                    }
                    rslen=0;
                }else{
                    rslen++;
                }
            }
        }
    }
}



int nrowBed(const char* fname, const char* reg){
    
    int i;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos1, pos2;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    int nrow=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d", chr, &pos1, &pos2);
        nrow++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}

int loadUKBB(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val){
    
    int i, k;
    int nrow;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    
    nrow = nrowBed(fname, reg);
    (*pchr)=(char*)calloc(1000, sizeof(char));
    (*ppos1)=(int*)calloc(nrow, sizeof(int));
    (*ppos2)=(int*)calloc(nrow, sizeof(int));
    (*val) = (double*)calloc(nrow, sizeof(double));
    char* al1; al1 = (char*)calloc(100, sizeof(char));
    char* al2; al2 = (char*)calloc(100, sizeof(char));
    char* rsid; rsid = (char*)calloc(100, sizeof(char));
    int nSmp;
    double ac, ytx, beta, se, tstat, pval;
    
    int gstart=270000000, gend=0;
    
    // chr, start, end, allele 1, allele 2, rsid, nCompleteSamples, AC, ytx, beta, se, tstat, pval
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    i=0;
    int nchar, ncharval;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", *pchr, (*ppos1)+i, (*ppos2)+i, al1, al2, rsid, &nSmp, &ac, &ytx, &beta, &se, &tstat, &pval);
        (*val)[i] = getLogWABFfromBetaSE(beta, se);
        i++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}

int loadBed1(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val, int nrow){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int isalloc=1;  // memory has been allocated
    if(nrow<0){
        isalloc=0;
        nrow = nrowBed(fname, reg);
        (*pchr)=(char*)calloc(1000, sizeof(char));
        (*ppos1)=(int*)calloc(nrow, sizeof(int));
        (*ppos2)=(int*)calloc(nrow, sizeof(int));
    }
    
    int gstart=270000000, gend=0;
    
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    i=0;
    int nchar, ncharval;
    int nval=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d%n", *pchr, (*ppos1)+i, (*ppos2)+i, &nchar);
        if(i==0){
            for(k=nchar; k<strlen(str.s); k++){
                if(str.s[k]=='\t'){nval++;}
            }
            if(isalloc==0)(*val) = (double*)calloc(nrow*nval, sizeof(double));
        }
        for(k=0; k<nval; k++){
            sscanf(str.s+nchar, "\t%lf%n", (*val)+i+k*nrow, &ncharval);
            nchar += ncharval;
        }
        i++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}

int loadBed(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val){
    return loadBed1(fname, reg, pchr, ppos1, ppos2, val, -1);
}

double getCovFromBedG(const char* fname, char* chr, int pos){
    char* reg; reg = (char*)calloc(1000, sizeof(char));
    sprintf(reg, "%s:%d-%d", chr, pos, pos+1);
    //fprintf(stderr, "%s %s\n", fname, reg);
    int n = nrowBed(fname, reg);
    int i;
    char* chrs; chrs = (char*)calloc(n, sizeof(char));
    int* pos1; pos1 = (int*)calloc(n, sizeof(int));
    int* pos2; pos2 = (int*)calloc(n, sizeof(int));
    double* val; val = (double*)calloc(n, sizeof(double));
    loadBed1(fname, reg, &chrs, &pos1, &pos2, &val, n);
    double res = val[0];
    free(reg);
    free(chrs);
    free(pos1);
    free(pos2);
    free(val);
    return res;
}

double cov2eta(double x, double beta0, double* beta1, double* xk, int nk, int EXPIT){
    double res = beta0 + beta1[0]*x;
    int i;
    //fprintf(stderr, "%lf %lf ", x, beta0);
    for(i=0; i<nk; i++){
        //fprintf(stderr, "%lf ", beta1[i+1]);
        res += beta1[i+1] * rk1(x, xk[i]);
    }
    //fprintf(stderr, "\n");
    return EXPIT>0 ? expit(res) : res;
}

int nrowVCF(const char* fname, const char* reg, int binarize, int* pnvars, int* psamplesize){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nal;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        
        if(nvars==0){nsamples = countFields(body, '\t');}
        
        nvars++;
        nal = countFields(a1, ',')+1;
        
        nbivars += achoose(nal);
    }
    tbx_itr_destroy(itr);
    (*pnvars) = nvars;
    (*psamplesize) = nsamples;
    if(binarize>0){
        return nbivars;
    }else{
        return nvars;
    }
}


int getBiVCF(const char* fname, const char* reg, double** pds1, int* psamplesize, int* pnbivars, char** pchr, int** ppos, char*** prss, char*** pba0, char*** pba1, int** vtype){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    (*pchr)   =(char*)calloc(1000, sizeof(char));
    char* chr; chr = (*pchr);
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info=(char*)calloc(1000, sizeof(char));
    VCF_info vinfo;
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    double* gl; gl = (double*)calloc(210, sizeof(double));
    double* ds; ds = (double*)calloc(20, sizeof(double));
    double* asc; asc = (double*)calloc(20, sizeof(double));
    double* ap; ap = (double*)calloc(40, sizeof(double));
    int*    gt; gt = (int*)calloc(2, sizeof(int));
    double* d;  d  = (double*)calloc(400, sizeof(double));
    
    int* formatID; formatID=(int*)calloc(100, sizeof(int));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0, nbivars1;
    int nal;
    int nfields;
    
    int ncol_ds1 = nrowVCF(fname, reg, 1, &nvars, &nsamples);
    //fprintf(stderr, "getBiVCF %d\n", ncol_ds1);
    (*pnbivars) = ncol_ds1;
    (*psamplesize) = nsamples;
    (*pds1) = (double*)calloc(nsamples*ncol_ds1, sizeof(double));
    (*pba0) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*pba1) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*prss) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*ppos) = (int*)calloc(ncol_ds1, sizeof(int));
    (*vtype) = (int*)calloc(ncol_ds1, sizeof(int));
    
    int* work; work=(int*)calloc(100, sizeof(int));
    char* cwork; cwork=(char*)calloc(1000, sizeof(char));
    int offs;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    int l=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        //if(l%100==0)fprintf(stderr, "%d\n", l);
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        body = str.s+nchar;
        nal = countFields(a1, ',')+1;
        
        nfields = parseFormat(format, formatID);
        parseInfo(info, cwork, &vinfo);
        
        nbivars1 = achoose(nal);
        //fprintf(stderr, "achoose=%ld\n", nbivars1);
        for(i=nbivars; i<nbivars+nbivars1; i++){
            (*pba0)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*pba1)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*prss)[i] = (char*)calloc(strlen(rs)+1, sizeof(char));
            strcpy((*prss)[i], rs);
            (*ppos)[i] = pos;
            (*vtype)[i] = vinfo.VT;
        }
        if(nal>2){
            getCombAlk(nal, 0, work, -1, nbivars1, a0, a1, (*pba0)+nbivars, (*pba1)+nbivars);
        }else{
            strcpy((*pba0)[nbivars], a0);
            strcpy((*pba1)[nbivars], a1);
        }
        for(i=nbivars; i<nbivars+nbivars1; i++){
            if((*vtype)[i]<0){
                if(strlen((*pba0)[i])==1 && strlen((*pba1)[i])==1){
                    (*vtype)[i] = VT_SNP;
                }else{
                    (*vtype)[i] = VT_INDEL;
                }
            }
        }
        offs = 0;
        //fprintf(stderr, "%s\n", body);
        for(i=0; i<nsamples; i++){
            //fprintf(stderr, "%s\n", body+offs);
            //int j; for(j=0; j<20; j++){fprintf(stderr, "%c", body[offs+j]);}fprintf(stderr, "\n");
            offs += parseCell(body+offs, nal, gt, gl, ap, asc, ds, formatID, nfields);
            //clearAs0(ds, nal);
            //ds[gt[0]]++;
            //ds[gt[1]]++;
            if(doseFormatExist(formatID, nfields, FORMAT_DS)==0){
                if(gt[0]<0){
                    ds[1]=-999999;
                }else{
                    if(doseFormatExist(formatID, nfields, FORMAT_GP)==1){
                        gl2ap(gt, gl, nal, ap, d);;
                    }else{
                        clearAs0(ap, nal*2);
                        ap[gt[0]]++;
                        ap[gt[1]+nal]++;
                    }
                    for(k=0; k<nal; k++){
                        ds[k] = ap[k]+ap[k+nal];
                    }
                }
            }
            if(nal>2){
                getCombk(nal, 0, work, -1, nbivars1, ds, (*pds1)+nbivars*nsamples+i, nsamples);
            }else{
                //getCombk(nal, 0, work, -1, nbivars1, ds, (*pds1)+nbivars*nsamples+i, nsamples);
                (*pds1)[nbivars*nsamples+i] = ds[1];
            }
        }
        nbivars += nbivars1;
        l++;
    }
    //fprintf(stderr, "finish\n");
    tbx_itr_destroy(itr);
    
    
    for(i=0; i<nbivars; i++){
        //fprintf(stderr, "%s %s %s\n", rss[i], ba0[i], ba1[i]);
    }
    
}

int getNpeaks(int* pcents, int n, int a, int b){
    int i;
    int res=0;
    for(i=0; i<n; i++){
        if(pcents[i]>a && pcents[i]<b){res++;}
    }
    return res;
}

// 0 .. npeaks in [-2W, 2W] for prior prob
// 0 .. M in [-W, W] for testing

int bayeslm(int argc, char** argv){
    
    int i, j, mode=MODE_N;
    
    // Verbose mode
    int verbose=0;
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--verbose")==0 || strcmp(argv[i], "-v")==0){verbose=1; break;}}
    
    // VCFs
    const char* fname = NULL;
    const char* fname2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--vcf")==0 || strcmp(argv[i], "-g")==0){fname = fname2 = argv[i+1]; break;}}
    if(fname==NULL){fprintf(stderr, "VCF file is missing.\n"); return 1;}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--vcf2")==0 || strcmp(argv[i], "-g2")==0){fname2  = argv[i+1]; break;}}
    if(verbose>0){fprintf(stderr, "VCF : %s %s\n", fname, fname2);};
    
    // cis-window
    int tss = 0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--window-centre")==0 || strcmp(argv[i], "-c")==0){tss = atoi(argv[i+1]); break;}}
    if(tss==0){fprintf(stderr, "Window centre is missing.\n"); return 1;}
    char* chrom=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--window-chromosome")==0 || strcmp(argv[i], "-s")==0){chrom = argv[i+1]; break;}}
    if(chrom==NULL){fprintf(stderr, "Chromosome is missing.\n"); return 1;}
    int wsize=500000;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--window-size")==0 || strcmp(argv[i], "-w")==0){wsize = atoi(argv[i+1])/2; break;}}
    
    // variant-level prior
    double* beta = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--variant-level")==0 || strcmp(argv[i], "-p1")==0){
        if(verbose>0){fprintf(stderr, "Beta1 : %s\n", argv[i+1]);};
        //int nbeta=gzfdscanf(argv[i+1], &beta);
        int nbeta=bdfscanf1h(argv[i+1], &beta, 10, 0);
        if(verbose>0)printV(beta, nbeta); 
        mode=MODE_P;
        break;}
    }
    double* beta2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--variant-level2")==0){
        if(verbose>0){fprintf(stderr, "Beta2 : %s\n", argv[i+1]);};  
        //int nbeta2=gzfdscanf(argv[i+1], &beta2); 
        int nbeta2=bdfscanf1h(argv[i+1], &beta2, 10, 0); 
        if(verbose>0)printV(beta2, nbeta2);
        break;}
    }
    if(beta2 == NULL){beta2 = beta;}
    
    // peak-pair-level prior
    double* beta_psi;
    double* Psi;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-pair-level")==0 || strcmp(argv[i], "-p3")==0){
        mode=MODE_PP;
        Psi=(double*)calloc(4, sizeof(double));
        if(verbose>0){fprintf(stderr, "Psi1 : %s\n", argv[i+1]);};
        bdfscanf1h(argv[i+1], &beta_psi, 12, 0);
        if(verbose>0){fprintf(stderr, "psi=");printV(beta_psi, 12);}
        break;}
    }
    
    // read counts
    const char* fnamey  = NULL;
    const char* fnamey2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--normalised-count")==0 || strcmp(argv[i], "-i")==0){fnamey = fnamey2 = argv[i+1]; break;}}
    if(fnamey==NULL){fprintf(stderr, "FPKM file is missing.\n"); return 1;}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--normalised-count2")==0 || strcmp(argv[i], "-i2")==0){fnamey2 = argv[i+1]; mode=MODE_C; break;}}
    if(verbose>0){fprintf(stderr, "FPKM : %s %s\n", fnamey, fnamey2);};
    
    
    
    // feature IDs
    int fid = 0;      // start point of fpkm, starting from 1
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-id")==0 || strcmp(argv[i], "-j")==0){fid = atoi(argv[i+1]); break;}}
    if(fid==0){fprintf(stderr, "Feature ID is missing.\n"); return 1;}
    int fid2 = fid+1;      // pair id, starting from 1
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-id2")==0 || strcmp(argv[i], "-j2")==0){fid2 = atoi(argv[i+1]); break;}}
    if(verbose>0){fprintf(stderr, "FID : %d %d\n", fid, fid2);};
    
    
    // GWAS
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--uk10k")==0){mode=MODE_G;}}
    
    // output file name
    gzFile outf=NULL; 
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--output")==0 || strcmp(argv[i], "-o")==0){outf = gzopen(argv[i+1], "ab6f");}}
    
    
    // 
    if(tss<1){tss=1;}
    
    
    // region in which variants are used
    char* reg; reg = (char*)calloc(1000, sizeof(char));
    int regstart, regend;
    if(mode==MODE_N){
        regstart = tss-wsize;
        regend   = tss+wsize;
    }else if(mode==MODE_P){
        regstart = tss-  wsize;
        regend   = tss+2*wsize;
    }else if(mode==MODE_C){
        regstart = tss-2*wsize;
        regend   = tss+2*wsize;
    }else if(mode==MODE_G){
        regstart = tss-wsize;
        regend   = tss+wsize;
    }else if(mode==MODE_PP){
        regstart = tss-2*wsize;
        regend   = tss+2*wsize;
    }
    if(regstart<=0){regstart = 1;}
    sprintf(reg, "%s:%d-%d", chrom, regstart, regend);
    if(verbose>0)fprintf(stderr, "Region=%s\n", reg);
    
    
    
    // flanking feature regions
    char* peakbed=NULL;
    int npeaks=0;
    int* pos1bed;
    int* pos2bed;
    int* pcents;
    //double* midp;
    int fid_bed=-1;
    int fid2_bed=0;
    int fid_bed_end=0;// end + 1
    int fid_bed_sta=0;
    double* ph; // relative peak height
    int M=0; // num of peaks for pairewise tests tss < peaks < tss+wsize are tested.
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-bed")==0 || strcmp(argv[i], "-p")==0){
        peakbed=argv[i+1]; 
        char* chrbed;
        npeaks = loadBed(peakbed, reg, &chrbed, &pos1bed, &pos2bed, &ph);
        pcents = (int*)calloc(npeaks, sizeof(int));
        //midp   = (double*)calloc(npeaks, sizeof(double));
        for(j=0; j<npeaks; j++){
            pcents[j] = (pos1bed[j]+pos2bed[j])/2;
            //midp[j] = ((double)pos1bed[j]+(double)pos2bed[j])/2.0;
            if(pcents[j]<tss-wsize){ fid_bed_sta=j+1; } // fprintf(stderr, "%d %d %d", fid_bed_sta, pcents[j], tss-wsize);}
            if(pcents[j]<tss+wsize){ fid_bed_end=j+1; }
            if(pos1bed[j]<=tss && tss<=pos2bed[j]){// target peak
                fid_bed = j;
            }else if(pos2bed[j]<tss){
                fid_bed = j;
            }
        }
        if(mode>0){
            if(mode==MODE_P){
                M = fid_bed_end - fid_bed - 1;
                fid2_bed=fid_bed+1;
            }else if(mode==MODE_C || mode==MODE_PP){
                M = fid_bed_end - fid_bed_sta;
                //fid2 -= (fid_bed-fid_bed_sta) + 1;
                fid2_bed = fid_bed_sta;
            }
            //if(verbose>0)fprintf(stderr, "fid_bed_sta=%d fid_bed_end=%d fid_bed=%d fid2_bed=%d fid2=%d", fid_bed_sta, fid_bed_end, fid_bed, fid2_bed, fid2);
        }
        break;
    }}
    //fprintf(stdout, "%d %d %d\n", pos1bed[fid2_bed], pos2bed[fid2_bed], fid2);
    if(verbose>0){fprintf(stderr, "mode=%d fid=%d fid2=%d fid_bed=%d fid2_bed=%d fid_bed_sta=%d, fid_bed_end=%d N peaks=%d\n", mode, fid, fid2, fid_bed, fid2_bed, fid_bed_sta, fid_bed_end, M);}
    if(M==0 && (mode==MODE_P || mode==MODE_C || mode==MODE_PP)){if(verbose>0){fprintf(stderr, "no paired features found\n");} if(outf!=NULL){gzclose(outf);}; return 0;}
    if(verbose>0){fprintf(stderr, "%d annotation features found.\n", npeaks);}
    
   
    
    
    
    // feature-level prior
    double* Pi1;    Pi1    = (double*)calloc(M+1, sizeof(double));
    double* Pi1p1;  Pi1p1  = Pi1 + 1;
    double* Pi1_a;  Pi1_a  = (double*)calloc(M+1, sizeof(double)); // not implemented for different trait 2
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-level")==0 || strcmp(argv[i], "-p2")==0){
        if(verbose>0){fprintf(stderr, "Feature-level prior : %s\n", argv[i+1]);}; 
        int nbeta=bdfscanf1(argv[i+1], Pi1, 1, fid-1);
        for(j=0; j<argc-1; j++){if(strcmp(argv[j], "--feature-level2")==0){// MODE_C
            if(verbose>0){fprintf(stderr, "Feature-level prior (trait 2) : %s\n", argv[j+1]);}; 
            nbeta=bdfscanf1(argv[j+1], Pi1+1, M, fid2-1);
            break;
        }}
        if(mode==MODE_P || mode==MODE_PP){ nbeta=bdfscanf1(argv[i+1], Pi1+1, M, fid2-1); }
        
        // Pi1_a for bothsides==0
        char* regtemp; regtemp = (char*)calloc(100, sizeof(char));
        sprintf(regtemp, "%s:%d-%d", chrom, tss-wsize, tss); // n of peaks on the 500Kb left hand side from peak j including j
        int np0 = getNpeaks(pcents, npeaks, tss-wsize, tss+2); //nrowBed(peakbed, regtemp);
        int np2;
        for(j=0; j<M; j++){
            //sprintf(regtemp, "%s:%d-%d", chrom, tss-wsize, (int)(midp[fid2_bed+j])+wsize); // n of peaks on the 500Kb right hand side from peak k including k
            //fprintf(stderr, "%d %s fid2=%d j=%d %lf \n", fid2_bed, regtemp, fid2_bed, j, midp[fid2_bed+j]+wsize);
            np2 = getNpeaks(pcents, npeaks, tss-wsize, pcents[fid2_bed+j]+wsize); //nrowBed(peakbed, regtemp);
            Pi1_a[j+1] = (bdfsum(argv[i+1], np2, fid-np0) - Pi1[0] - Pi1[j+1])/((double)(np2-2));
            if(isnan(Pi1_a[j+1])>0 || Pi1_a[j+1]<0.0 || Pi1_a[j+1]>1.0){Pi1_a[j+1] = 1e-20;}
            //fprintf(stderr, "%d %d %d \n", j+1, np2, fid-np0);
            //Pi1_a = nk_mean(Pi1+1, M);
            //return 0;
        }
        break;
    }}
    if(verbose>0){
        fprintf(stderr, "Pi1:  "); printV2(Pi1, M+1);
        fprintf(stderr, "Pi1_a:"); printV2(Pi1_a,M+1);
    }
    
    
    
    
    // loading genotypes from vcf1
    double* ds;
    int nbivars, samplesize;
    char* chr;
    int* pos;
    char** ba0;
    char** ba1;
    char** rss;
    int* vt;
    if(verbose>0){fprintf(stderr, "Loading genotype dose from %s in %s.\n", fname, reg);}
    getBiVCF(fname, reg, &ds, &samplesize, &nbivars, &chr, &pos, &rss, &ba0, &ba1, &vt);
    if(verbose>0){fprintf(stderr, "Genotype dose has been successfully loaded from %s in %s (%d x %d).\n", fname, reg, nbivars, samplesize);}
    
    double maf0 = 0.05;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--minor-allele-frequency")==0){maf0 = (double)atof(argv[i+1]); break;}}
    
    double* af;     af     = (double*)calloc(nbivars, sizeof(double));
    double* rsq;    rsq    = (double*)calloc(nbivars, sizeof(double));
    double* bf;     bf     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* bfmr;   bfmr   = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* bfmr2;  bfmr2  = (double*)calloc(nbivars * (M+1), sizeof(double));
    
    // feature-level prior was here
    int*    loccat; loccat = (int*)calloc(nbivars * (M+1), sizeof(int));
    int*    loccatid; loccatid = (int*)calloc(nbivars, sizeof(int)); for(i=0;i<nbivars;i++){loccatid[i]=-1;}
    double* w;      w      = (double*)calloc(nbivars, sizeof(double)); // 1.0: effective loci; 0.0: unused (MAF==0 etc.)
    
    // vcf2
    double* ds2;
    int nbivars2, samplesize2;
    char* chr2;
    int* pos2;
    char** ba02;
    char** ba12;
    char** rss2;
    int* vt2;
    double* af2;
    double* rsq2;
    
    if(mode==MODE_C){
        if(verbose>0){fprintf(stderr, "Loading genotype dose from %s in %s.\n", fname2, reg);}
        getBiVCF(fname2, reg, &ds2, &samplesize2, &nbivars2, &chr2, &pos2, &rss2, &ba02, &ba12, &vt2);
        if(verbose>0){fprintf(stderr, "Genotype dose has been successfully loaded from %s in %s (%d x %d).\n", fname2, reg, nbivars2, samplesize2);}
        af2     = (double*)calloc(nbivars, sizeof(double));
        rsq2    = (double*)calloc(nbivars, sizeof(double));
        if(nbivars != nbivars2){fprintf(stderr, "Different VCF conposition!\n"); if(outf!=NULL){gzclose(outf);}; return 1;}
    }else{
        ds2 = ds;
        samplesize2 = samplesize;
        nbivars2 = nbivars;
        chr2 = chr;
        pos2 = pos;
        ba02 = ba0;
        ba12 = ba1;
        rss2 = rss;
        vt2  = vt;
        af2  = af;
        rsq2 = rsq;
    }
    if(verbose>0){fprintf(stderr, "\n");}
    
    
    
    
    
    // posterior prob from pairwise model
    double* Zj;  Zj  = (double*)calloc(nbivars*2, sizeof(double));
    
    
    
    
    
    // loading read counts
    FILE* fy;
    fy = fopen(fnamey, "rb");
    double* y;  y  = (double*)calloc(samplesize, sizeof(double));
    fseek(fy, samplesize*(fid-1)*sizeof(double), SEEK_SET);
    int fread_info = fread(y, sizeof(double), samplesize, fy);
    // loading read counts for paired trait
    double* y2;
    if(verbose>0){fprintf(stderr, "%d feature fpkms are loaded from %s.\n", samplesize, fnamey);}
    if(mode==MODE_P || mode==MODE_C || mode==MODE_PP){
        FILE* fy2;
        fy2 = fopen(fnamey2, "rb");
        y2 = (double*)calloc(samplesize2 * M, sizeof(double));
        fseek(fy2, samplesize2*(fid2-1)*sizeof(double), SEEK_SET);
        fread_info = fread(y2, sizeof(double), samplesize2 * M, fy2);
        if(verbose>0){fprintf(stderr, "%d x %d feature fpkms are loaded from %s.\n", M, samplesize2, fnamey2);}
    }
    if(verbose>0){fprintf(stderr, "\n");}
    
    // randomise y and y2
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--random-permutation")==0 || strcmp(argv[i], "-r")==0){
        srand((unsigned)(time(NULL)+getpid()));
        if(mode==MODE_N || mode==MODE_G){
            randomise(y, samplesize);
        }else if(mode==MODE_P || mode==MODE_C || mode==MODE_PP){
            randomise2(y, y2, M, samplesize);
        }
        break;
    }}
   
    // GxE
    double* workForBFCalc; workForBFCalc = (double*)calloc(10*samplesize + (3*4+2)*4 + 1, sizeof(double));
    double* env = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--env")==0){
        //fprintf(stderr, "Interaction: %s\n", argv[i+1]);
        FILE* fenv; fenv = fopen(argv[i+1], "rb");
        env = (double*)calloc(samplesize, sizeof(double));
        fread_info = fread(env, sizeof(double), samplesize, fenv);
        fclose(fenv);
        break;
    }}
 
    // computation of BFs
    double* be;     be     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* se;     se     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* eta;    eta    = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* eta0;   eta0   = (double*)calloc(nbivars, sizeof(double));
    double* covs;   covs   = (double*)calloc(nbivars*2, sizeof(double));
    int k, l, maxid;
    double maxbf = -1.0e10;
    int ntested=0;
    char* fcoverage=NULL; for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--coverage")==0){fcoverage = argv[i+1]; if(verbose>0){fprintf(stderr, "Coverage file: %s\n", fcoverage);}; break;}}
    if(verbose>0) fprintf(stderr, "Prior calculation...");
    for(j=0; j<nbivars; j++){
        af[j]  = nk_dsum(ds+j*samplesize, samplesize, 1)/2.0/(double)samplesize;
        rsq[j] = nk_var(ds+j*samplesize, ds+j*samplesize, samplesize)/af[j]/(1.0-af[j])/2.0;
        if(beta!=NULL){
            eta[j] = eta0[j] = (vt[j]>0 ? beta[0] : 0.0); 
            for(l=0; l<M; l++){ eta[j+(l+1)*nbivars] = (vt[j]>0 ? beta2[0] : 0.0);} // beta2[vt[j]]; }
        }
        //fprintf(stderr, "%lf %lf\n", af[j], rsq[j]);
        if(af[j]>maf0 && af[j]<1.0-maf0 && rsq[j]>0.3){
            //fprintf(stderr, ".");
            //Prior calculation
            for(k=0; k<npeaks; k++){// annnotation
                if(pos1bed[k]<=pos[j] && pos[j]<=pos2bed[k]){// in the kth peak
                    //fprintf(stderr,"+");
                    //double covAtj;
                    if(fcoverage!=NULL){ covs[j] = getCovFromBedG(fcoverage, chrom, pos[j]); } //fprintf(stderr, "%lf\n", covs[j]); }
                    if(beta!=NULL) eta0[j] += beta[2]; //cov2eta(covs[j]/ph[k], beta[4], beta+5, xk, 4, 0);
                    if(k==fid_bed){// is target peak?
                        loccat[j] = 1; loccatid[j]=k;
                        if(beta!=NULL) eta[j] += beta[1]; //cov2eta(covs[j]/ph[k], beta[3], beta+5, xk, 4, 0);
                    }else{
                        loccat[j] = 2; loccatid[j]=k;
                        if(beta!=NULL) eta[j] += beta[2]; //cov2eta(covs[j]/ph[k], beta[4], beta+5, xk, 4, 0);
                    }
                    for(l=0; l<M; l++){// pairwise
                        if(k==l+fid2_bed){
                            loccat[j+(l+1)*nbivars] = 1;
                            if(beta2!=NULL) eta[j+(l+1)*nbivars] += beta2[1]; //cov2eta(covs[j]/ph[k], beta2[3], beta2+5, xk, 4, 0);
                        }else{
                            loccat[j+(l+1)*nbivars] = 2;
                            if(beta2!=NULL) eta[j+(l+1)*nbivars] += beta2[2]; //cov2eta(covs[j]/ph[k], beta2[4], beta2+5, xk, 4, 0);
                        }
                    }
                    break;
                }
            }
            
            // BF calculation
            //fprintf(stderr, "%s\t%lf\n", rss[j],ds[j*samplesize]);
            w[j] = 1.0;
            //bf[j] = getLogBF(ds+j*samplesize, y, samplesize, sqrt(10.0), work);
            if(env==NULL){
                bf[j] = getLogWABF(ds+j*samplesize, y, samplesize);
                be[j] = getBeta(ds+j*samplesize, y, samplesize, se+j);
            }else{
                bf[j] = getLogWABFInter(ds+j*samplesize, y, env, samplesize, workForBFCalc);
            }
            if(mode==MODE_P || mode==MODE_C || mode==MODE_PP){// for pairwise
                af2[j]  = nk_dsum(ds2+j*samplesize2,       samplesize2, 1)/2.0/(double)samplesize2;
                rsq2[j] = nk_var( ds2+j*samplesize2, ds2+j*samplesize2, samplesize2)/af2[j]/(1.0-af2[j])/2.0;
                if(af2[j]>maf0 && af2[j]<1.0-maf0 && rsq2[j]>0.3 && af[j]>maf0 && af[j]<1.0-maf0 && rsq[j]>0.3){
                    for(k=0; k<M; k++){
                        //bf[j+(k+1)*nbivars] = getLogBF(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, 10.0, work);
                        bf[j+(k+1)*nbivars]    = getLogWABF(  ds2+j*samplesize2,    y2+k*samplesize2, samplesize2);
                        if(ds==ds2 || samplesize==samplesize2){
                            bfmr[j+(k+1)*nbivars]  = getLogWABFMR(ds2+j*samplesize2, y, y2+k*samplesize2, samplesize2);
                            bfmr2[j+(k+1)*nbivars] = getLogWABFMR(ds2+j*samplesize2, y2+k*samplesize2, y, samplesize2);
                        }else{
                            // fprintf(stderr, "atac-eqtl");
                            // expression - atac : not implemented
                            bfmr[ j+(k+1)*nbivars] = getLogWABFMR2(ds+ j*samplesize,  y,                samplesize,  ds2+j*samplesize2, y2+k*samplesize2, samplesize2);
                            bfmr2[j+(k+1)*nbivars] = getLogWABFMR2(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, ds+j*samplesize,   y,                samplesize );
                        }
                        be[j+(k+1)*nbivars]    = getBeta(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, se+(j+(k+1)*nbivars));
                    }
                }else{
                    w[j] = 0.0;
                }
            }
            //if(strcmp("rs2409780",rss[j])==0){w[j]=0.0;}
            //if(loccatid[j]==106){fprintf(stderr, "%s %d %d %lf\n", rss[j], loccat[j], loccatid[j], w[j]);}
            if(maxbf < bf[j]){maxid = j; maxbf=bf[j];}
            ntested++;
            if(mode==MODE_N){
                if(outf==NULL){
                    printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]);
                    for(k=0; k<M; k++){
                        printf("%d\t",  loccat[j+k*nbivars]);
                        printf("%lf\t%lf\t%lf\t%lf\t", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars], eta[j+k*nbivars]);
                    }
                    printf("%d\t", loccat[j+k*nbivars]);
                    printf("%lf\t%lf\t%lf\t%lf\t%lf\n", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars], eta[j+k*nbivars], covs[j]);
                }else{
                    gzprintf(outf, "%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]);
                    for(k=0; k<M; k++){
                        gzprintf(outf, "%d\t",  loccat[j+k*nbivars]);
                        gzprintf(outf, "%lf\t%lf\t%lf\t", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars]);
                    }
                    gzprintf(outf, "%d\t", loccat[j+k*nbivars]);
                    gzprintf(outf, "%lf\t%lf\t%lf\n", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars]);
                }
            }
        }
    }
    if(verbose>0) fprintf(stderr, "Done\n");
    
    // pairwise hierahical model for GWAS
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--uk10k")==0){
        FILE* fbf2list; fbf2list = fopen(argv[i+1], "r");
        char* fbf2; fbf2  = (char*)calloc(1000, sizeof(char));
        
        while(fscanf(fbf2list, "%s\n", fbf2)!=EOF){
            //fprintf(stderr, "%s\n", fbf2);
        
            int nrowbf2 = nrowBed(fbf2, reg);
            char* pfbf2; for(k=strlen(fbf2)-1; k>=0; k--){if(fbf2[k]=='/'){pfbf2 = fbf2+k+1; break;}}
                //if(nrowbf2==0){fprintf(stderr, "no gwas bfs\n"); return 0;}
            if(nrowbf2<10){
                fprintf(stderr, "No line in %s in %s\n", fbf2, reg);
            }else{
                char* chrbf2;    chrbf2  = (char*)calloc(1000, sizeof(char));
                int* posbf2;     posbf2  = (int*)calloc(nrowbf2, sizeof(int));
                double* bf2orig; bf2orig = (double*)calloc(nrowbf2, sizeof(double));
                double* bf2;     bf2     = (double*)calloc(nbivars, sizeof(double));
                int* pos2bf2; pos2bf2 = (int*)calloc(nrowbf2, sizeof(int));
                loadBed(fbf2, reg, &chrbf2, &posbf2, &pos2bf2, &bf2orig);
                //loadUKBB(fbf2, reg, &chrbf2, &posbf2, &pos2bf2, &bf2orig);
                
                //printV(bf2orig, nrowbf2);
                clearAs0(w, nbivars);
                expandInt(pos2bf2, bf2orig, nrowbf2, pos, bf2, nbivars, w);
                //fprintf(stderr, "N ol vars=%lf\n", nk_dsum(w,nbivars,1)); 
                
                //double* pp3;  pp3  = (double*)calloc(3,  sizeof(double));
                double* pp12; pp12 = (double*)calloc(12, sizeof(double));
                //pwhm13(bf, bf2, bf, vt, loccat, w, nbivars, beta, pp3, pp12, &phi0, &del0);
                // Param 5
                pwhmNewAtacGwas(bf, bf2, Pi1[0], w, nbivars, pp12);
                //printf("%d", fid);
                //for(j=0; j<6; j++){printf("\t%lf", log(pp12[j]));}printf("\n");
                
                
                //gzFile outf=NULL;
                //for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
                if(outf==NULL){
                    printf("%d\t%s\t%d", nrowbf2, pfbf2, fid);
                    for(k=0; k<6; k++){printf("\t%lf", log(pp12[k]));} printf("\n"); 
                }else{
                    gzprintf(outf, "%d\t%s\t%d", nrowbf2, pfbf2, fid);
                    for(k=0; k<6; k++){gzprintf(outf, "\t%lf", log(pp12[k]));} gzprintf(outf, "\n"); 
                }
            }
        }
        if(outf!=NULL){gzclose(outf);}
        return 0;
    }}
    
   
    // pairewise hierahical model for atac & eqtl
    if(mode==MODE_C || mode==MODE_P || mode==MODE_PP){
        if(verbose>0)fprintf(stderr, "Pairwise model\n");
        double* pp13; pp13 = (double*)calloc(12, sizeof(double));
        for(j=0; j<M; j++){
            int nloci = 0; // #vars <- 500K - P1 - p2 - 500K ->
            int geta = 0;
            int cis_sta = mini(tss, pcents[fid2_bed+j]) - wsize; if(cis_sta<0){cis_sta=1;}
            int cis_end = maxi(tss, pcents[fid2_bed+j]) + wsize;
            for(k=0; k<nbivars; k++){
                if(pos[k]  <  cis_sta){geta++;}
                if(cis_sta <= pos[k] && pos[k] <= cis_end){nloci++;}
                if(cis_end <  pos[k]){break;}
            }
            clearAs0(pp13, 13);
            
            if(mode==MODE_P || mode==MODE_PP){ // same trait
                if(fid!=fid2+j){
                    if(mode==MODE_PP){// posterior calculation after fitting
                        if(verbose>0){fprintf(stderr, "Posterior prob calculation and lead variant detection\n");}
                        
                        double xk[4] = {0.1470588, 0.3823529, 0.6176471, 0.8529412};// for spline
                        
                        // prior Psi
                        Psi[0] = 0.0;
                        Psi[1] = beta_psi[0];
                        Psi[2] = beta_psi[6];
                        Psi[3] = beta_psi[6];
                        
                        double apeakdis = ((double)abs(pcents[fid_bed]-pcents[fid2_bed+j]))/500000.0;
                        
                        Psi[1] += apeakdis*beta_psi[1];
                        Psi[2] += apeakdis*beta_psi[7];
                        Psi[3] += apeakdis*beta_psi[7];
                        
                        for(k=0; k<4; k++){
                            Psi[1] += rk1(apeakdis, xk[k])*beta_psi[2+k];
                            Psi[2] += rk1(apeakdis, xk[k])*beta_psi[8+k];
                            Psi[3] += rk1(apeakdis, xk[k])*beta_psi[8+k];
                        }
                        double tot_psi = exp(Psi[0]) + exp(Psi[1]) + exp(Psi[2]) + exp(Psi[3]);
                        for(k=0; k<4; k++){ Psi[k] = exp(Psi[k])/tot_psi; }
                        
                        //fprintf(stderr, "%d %d %d %d, %lf, %lf %lf %lf\n", fid, fid2+j, pcents[fid_bed], pcents[fid2_bed+j], Psi[0], Psi[1], Psi[2], Psi[3]);
                        
                        // variant location relative to peak regions
                        geta = 0;
                        nloci = nbivars;
                        pwhmfm( bf+geta, bf+(j+1)*nbivars+geta, bfmr+(j+1)*nbivars+geta, bfmr2+(j+1)*nbivars+geta, eta0+geta, eta+geta, eta+(j+1)*nbivars+geta, Pi1[0], Pi1[j+1], w+geta, nloci, Psi, Zj+geta);
                    }else if(mode==MODE_P){
                        if(verbose>0){fprintf(stderr, "Pi1=%lf Pi1_pair=%lf Pi1_a=%lf\n", Pi1[0], Pi1[j+1], Pi1_a[j+1]);}

                        
                        pwhmnew(bf+geta, bf+(j+1)*nbivars+geta, bfmr+(j+1)*nbivars+geta, bfmr2+(j+1)*nbivars+geta, eta0+geta, eta+geta, eta+(j+1)*nbivars+geta, Pi1[0], Pi1[j+1], Pi1_a[j+1], w+geta, nloci, pp13, loccatid+geta);
                            
                        double apeakdis = ((double)abs(pcents[fid_bed]-pcents[fid2_bed+j]))/500000.0;
                        
                        if(outf==NULL){
                            printf("%d\t%d\t%lf\t", fid, fid2+j, apeakdis);
                            for(k=0; k<9; k++){printf("%lf ", log(pp13[k]));} printf("%lf\n", log(pp13[k])); 
                        }else{
                            gzprintf(outf, "%d\t%d\t%lf\t", fid, fid2+j, apeakdis);
                            for(k=0; k<9; k++){gzprintf(outf, "%lf ", log(pp13[k]));} gzprintf(outf, "%lf\n", log(pp13[k])); 
                        }
                    
                    }
                    
                }
            }else if(mode==MODE_C){ // different traits (e.g., atac-eqtl)
                coloc(bf+geta, bf+(j+1)*nbivars+geta, eta0+geta, Pi1[0], Pi1[j+1], w+geta, nloci, pp13);

                if(outf==NULL){
                    printf("%d\t%d\t", fid, fid2+j);
                    for(k=0; k<5; k++){printf("%lf ", log(pp13[k]));}printf("%lf\n", log(pp13[k]));
                }else{
                    gzprintf(outf, "%d\t%d\t", fid, fid2+j);
                    for(k=0; k<5; k++){gzprintf(outf, "%lf ", log(pp13[k]));}
                    gzprintf(outf, "%lf\n", log(pp13[5]));
                }
            }
        }
        if(mode==MODE_PP){ // pwhmfm() posterior probability estimate using causal inference information
            int cis_sta = tss - wsize; if(cis_sta<0){cis_sta=1;}
            int cis_end = tss + wsize;
            double totzj=0.0, totzjnom=0.0;
            double* varloc; varloc=(double*)calloc(3,sizeof(double));
            for(i=0; i<nbivars; i++){
                if(cis_sta <= pos[i] && pos[i] <= cis_end){varloc[loccat[i]] += Zj[i]/totzj;}
            }
            
            // no pairwise solution
            pwhmfm0(bf, eta, Pi1[0], w, nbivars, Zj+nbivars);
            for(i=0; i<nbivars; i++){
                if(cis_sta <= pos[i] && pos[i] <= cis_end){
                    totzj += Zj[i]; 
                    totzjnom += Zj[i+nbivars];
                }
            }

            //gzFile outf=NULL;
            //for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
            if(outf==NULL){
                if(1==1){
                    for(k=0; k<nbivars; k++){
                        if(af[k]>maf0 && af[k]<1.0-maf0 && rsq[k]>0.3){
                            printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n", fid, chr, pos[k], rss[k], ba0[k], ba1[k], af[k], rsq[k], vt[k], loccat[k], be[k], se[k], log(Zj[k]/totzj), log(Zj[k+nbivars]/totzjnom));
                        }
                    }
                }else{
                    printf("%d\t%lf\t%lf\t%lf\n", fid, varloc[0], varloc[1], varloc[2]);
                }
            }else{
                if(1==1){
                    for(k=0; k<nbivars; k++){
                        if(af[k]>maf0 && af[k]<1.0-maf0 && rsq[k]>0.3){
                            gzprintf(outf, "%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n", fid, chr, pos[k], rss[k], ba0[k], ba1[k], af[k], rsq[k], vt[k], loccat[k], be[k], se[k], log(Zj[k]/totzj), log(Zj[k+nbivars]/totzjnom));
                        }
                    }
                    //gzclose(outf);
                }else{
                    gzprintf(outf, "%d\t%lf\t%lf\t%lf\n", fid, varloc[0], varloc[1], varloc[2]);
                    //for(i=0; i<nbivars; i++){ if(cis_sta <= pos[i] && pos[i] <= cis_end){gzprintf(outf, "%lf\n", Zj[i]/totzj);} }
                    //gzclose(outf);
                }
            }
        }
    }
    if(outf!=NULL){ gzclose(outf);}
    return 0;
}







int main(int argc, char** argv){
    if(argc==1){usage_bayeslm(); return 1;}
    exp_gt_gtdsgl = 0;
    if(verbose_loadVCF>0)fprintf(stderr, "\n\nbayesLm\n\n");
    bayeslm(argc, argv);

}
