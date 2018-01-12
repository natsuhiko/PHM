#include <stdio.h>
#include <math.h> 
#include <gsl/gsl_cdf.h>


double getLogBF(double* g, double* y, int n, double sigma, double* work);

void getLogBFForR(double* g, double* y, int* n, double* sigma, double* work, double* lbf){
    (*lbf) = getLogBF(g, y, *n, *sigma, work);
}

double sign(double x){
    return x>0 ? 1.0 : -1.0; 
}

double ptqn(double tstat, double n){
    //return tstat;
    return sign(tstat) * gsl_cdf_ugaussian_Pinv( gsl_cdf_tdist_P(-fabs(tstat), n-2.) );
}

double getBeta(double* g, double* y, int n, double* se);

double getBeta(double* g, double* y, int n, double* se){
    int i;
    double mg, vg, my, vgy, vy, rgy2;
    mg = vg = my = vgy = vy = 0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){
        mg += g[i];
        my += y[i];
        vg += g[i]*g[i];
        vy += y[i]*y[i];
        vgy += g[i]*y[i];
    }
    mg /= dn;
    my /= dn;
    vg = vg/dn - mg*mg;
    vy = vy/dn - my*my;
    vgy= vgy/dn - my*mg;
    rgy2= vgy*vgy/vg/vy;
    
    (*se) = sqrt(vy/vg*(1.-rgy2)/(dn-2.));
    return vgy/vg;
}


double getLogWABF(double* g, double* x, int n){
    int i, j, k;
    double dn = (double)n;
    double g1, g2, x1, x2, gx, rgx2;
    double s2 = 10.;
    g1=g2=x1=x2=gx=rgx2=0.0;
    for(i=0; i<n; i++){
        x1 += x[i];
        gx += x[i]*g[i];
        x2 += x[i]*x[i];
        g1 += g[i];
        g2 += g[i]*g[i];
    }
    x1 /= dn;
    x2 /= dn;
    g1 /= dn;
    g2 /= dn;
    gx /= dn;
    x2 -= x1*x1;
    g2 -= g1*g1;
    gx -= g1*x1;
    rgx2 = gx*gx/g2/x2;
    
    //return x2/g2*(1.-rgx2)/(dn-2.);
    double lbf = 0.0;
    double r;
    double z2 = ptqn(sqrt((dn-2.)*rgx2/(1.-rgx2)), dn); z2=z2*z2;
    //double z2 = (dn-2.)*rgx2/(1.-rgx2);
    for(i=0; i<1; i++){
        r  = s2 / (s2 + x2/g2*(1.-rgx2)/(dn-2.));
        lbf += 0.5*log(1.-r) + 0.5*z2*r;
        s2 *= 2.0;
    }
    return lbf;
}

double getLogWABFfromBetaSE(double beta, double se){
    int i, j, k;
    double s2[3] = {0.01, 0.1, 0.5};
    
    double bf = 0.0;
    double r;
    double z = beta/se;
    for(i=0; i<3; i++){
        r   = s2[i] / (s2[i] + se*se);
        bf += exp( 0.5*log(1.-r) + 0.5*z*z*r ) / 3.0;
    }
    return log(bf);
}

double getLogWABFMR(double* g, double* x, double* y, int n){
    int i, j, k;
    double dn = (double)n;
    double g1, g2, x1, x2, y1, y2, gx, gy, xy;
    double s2 = 10.;
    g1=g2=x1=x2=y1=y2=gy=xy=gx=0.0;
    for(i=0; i<n; i++){
        x1 += x[i];
        x2 += x[i]*x[i];
        y1 += y[i];
        y2 += y[i]*y[i];
        g1 += g[i];
        g2 += g[i]*g[i];
        
        gx += x[i]*g[i];
        xy += x[i]*y[i];
        gy += g[i]*y[i];
    }
    x1 /= dn;
    x2 /= dn;
    y1 /= dn;
    y2 /= dn;
    g1 /= dn;
    g2 /= dn;
    gx /= dn;
    gy /= dn;
    xy /= dn;
    x2 -= x1*x1;
    y2 -= y1*y1;
    g2 -= g1*g1;
    
    gx -= g1*x1;
    gy -= g1*y1;
    xy -= x1*y1;
    
    gx /= sqrt(g2*x2);
    gy /= sqrt(g2*y2);
    xy /= sqrt(x2*y2);
    
    //fprintf(stderr, "%lf %lf %lf ", gx, gy, xy);
    
    double lbf = 0.0;
    double r;
    double z2 = ptqn(sqrt((dn-2.)*gy*gy/(1.-xy*xy+(gy/gx-xy)*(gy/gx-xy))), dn); z2 = z2*z2;
    for(i=0; i<1; i++){
        r = (1.-xy*xy+(gy/gx-xy)*(gy/gx-xy)) / ((dn-2.)*gx*gx) * y2/x2;
        r = s2 / (s2+r);
        lbf += 0.5*log(1.-r) + 0.5*z2*r;
        s2 *= 2.0;
    }
    //if(isnan(lbf)>0){fprintf(stderr, "%lf %lf %lf ", lbf, r, z2);}
    return lbf;
}


double getLogWABFMR2(double* g1, double* x, int n1, double* g2, double* y, int n2){
    int i, j, k;
    double dn1 = (double)n1;
    double dn2 = (double)n2;

    double g11, g12, g21, g22, x1, x2, y1, y2, gx, gy, z2gx, z2gy, beta2;
    double s2 = 100.;
    g11=g12=g21=g22=x1=x2=y1=y2=gy=gx=0.0;
    for(i=0; i<n1; i++){
        x1  += x[i];
        x2  += x[i]*x[i];
        g11 += g1[i];
        g12 += g1[i]*g1[i];
        gx  += g1[i]*x[i];
    }
    for(i=0; i<n2; i++){
        y1  += y[i];
        y2  += y[i]*y[i];
        g21 += g2[i];
        g22 += g2[i]*g2[i];
        gy  += g2[i]*y[i];
    }
    
    x1  /= dn1;
    x2  /= dn1;
    g11 /= dn1;
    g12 /= dn1;
    gx  /= dn1;
    x2  -= x1*x1;
    g12 -= g11*g11;
    gx  -= g11*x1;
    
    y1  /= dn2;
    y2  /= dn2;
    g21 /= dn2;
    g22 /= dn2;
    gy  /= dn2;
    y2  -= y1*y1;
    g22 -= g21*g21;
    gy  -= g21*y1;
    
    beta2 = gy*gy/g22/g22/gx/gx*g12*g12;
    
    gx  /= sqrt(g12*x2);
    z2gx = (dn1-2.) * gx*gx/(1.-gx*gx);
    gy  /= sqrt(g22*y2);
    z2gy = (dn2-2.) * gy*gy/(1.-gy*gy);
    
    //fprintf(stderr, "%lf %lf %lf %lf\n", z2gx, z2gy, x2, y2);
    //fprintf(stderr, "%lf %lf %lf ", gx, gy, xy);
    
    double lbf = 0.0;
    double r0  = beta2*(1./z2gx + 1./z2gy);
    double z2 = z2gx * z2gy / (z2gx + z2gy);
    double r;
    
    for(i=0; i<1; i++){
        r = s2 / (s2+r0);
        lbf += 0.5*log(1.-r) + 0.5*z2*r;
        s2 *= 2.0;
    }
    //fprintf(stderr, "%lf %lf %lf %lf\n", z2, r0, r, lbf);
    //if(isnan(lbf)>0){fprintf(stderr, "%lf %lf %lf ", lbf, r, z2);}
    return lbf;
}

// work length 6
double getLogBF(double* g, double* y, int n, double sigma, double* work){
    int i, j, k;
    double* invOmega; invOmega=work;
    double dn = (double)n;
    invOmega[0] = dn;
    invOmega[1] = 0.0;
    invOmega[2] = 1.0/sigma/sigma;
    double* xy; xy = work+3;
    xy[0] = xy[1] = xy[2] = 0.0;
    for(i=0; i<n; i++){
        invOmega[1] += g[i];
        invOmega[2] += g[i]*g[i];
        xy[0] += y[i];
        xy[1] += y[i]*g[i];
        xy[2] += y[i]*y[i];
    }
    double det = invOmega[0]*invOmega[2] - invOmega[1]*invOmega[1];
    double tmp  =  invOmega[2]/det;
    invOmega[2] =  invOmega[0]/det;
    invOmega[1] = -invOmega[1]/det;
    invOmega[0] =  tmp;
    
    double bxb;
    for(j=0; j<2; j++){
        for(k=0; k<2; k++){
            bxb += xy[j]*xy[k]*invOmega[j+k];
        }
    }
    
    return -0.5*log(det) + 0.5*log(dn) - log(sigma) - (dn/2.0)*(log(xy[2] - bxb) - log(xy[2] - xy[0]*xy[0]/dn) );
}

void getLogWABFforR(double* g, double* x, int* pn, double* bf){
    (*bf) = getLogWABF(g, x, *pn);
}
void getLogWABFMRforR(double* g, double* x, double* y, int* pn, double* bf){
    (*bf) = getLogWABFMR(g, x, y, *pn);
}
void getLogWABFMR2forR(double* g1, double* x, int* pn1, double* g2, double* y, int* pn2, double* bf){
    (*bf) = getLogWABFMR2(g1, x, *pn1, g2, y, *pn2);
}







