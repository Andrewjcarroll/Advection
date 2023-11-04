#ifndef KIM_FILTER_H
#define KIM_FILTER_H
#include <cmath>

double kc = 0.88;
double eps = 0.25;

#define IDX(i,j) ( (i) + n * (j) )

// kim_filter_cal_coeff {{{
void kim_filter_cal_coeff(double *c, double kc)
{
    double AF = 30.0 - 5.0*cos(kc) + 10*cos(2.0*kc) - 3.0*cos(3.0*kc);
    double alphaF = -(30.0*cos(kc) + 2.0*cos(3.0*kc)) / AF;
    double betaF = (18.0 + 9.0*cos(kc) + 6.0*cos(2.0*kc) - cos(3.0*kc)) 
                         / (2.0*AF);
    c[0] = AF;
    c[1] = alphaF;
    c[2] = betaF;
}
// }}}

// kim_filter_PQ {{{
void kim_filter_PQ(double *P, double *Q, double kc, double eps, int n)
{
    double c0[3];
    double cd[3];
    double cdd[3];
    double cddd[3];

    double t2 = sin(M_PI/2.0);
    double t3 = sin(M_PI/3.0);
    double t6 = sin(M_PI/6.0);

    double kcd = kc*(1.0 - eps*t6*t6);
    double kcdd = kc*(1.0 - eps*t3*t3);
    double kcddd = kc*(1.0 - eps*t2*t2);

    kim_filter_cal_coeff(c0, kc);
    kim_filter_cal_coeff(cd, kcd);
    kim_filter_cal_coeff(cdd, kcdd);
    kim_filter_cal_coeff(cddd, kcddd);

    const double t1 = cos(0.5*kc);
    const double aF1 = 30.0*t1*t1*t1*t1/c0[0];
    const double aF2 = -2.0*aF1/5.0;
    const double aF3 = aF1/15.0;
    const double aF0 = -2.0*(aF1 + aF2 + aF3);

    const double alphaF = c0[1];
    const double betaF  = c0[2];

    const double alphaFd = cd[1];
    const double betaFd  = cd[2];
    const double alphaFdd = cdd[1];
    const double betaFdd  = cdd[2];
    const double alphaFddd = cddd[1];
    const double betaFddd  = cddd[2];

    const double t1d = cos(0.5*kcd);
    const double aF1d = 30.0*t1d*t1d*t1d*t1d/cd[0];
    const double aF2d = -2.0*aF1d/5.0;
    const double aF3d = aF1d/15.0;
 
    const double BF = (1.0 - betaFdd)*(1.0 + 6.0*betaFdd + 60.0*betaFdd*betaFdd) 
                    + (5.0 + 35*betaFdd -29.0*betaFdd*betaFdd)*alphaFdd 
                    + (9.0 - 5.0*betaFdd)*alphaFdd*alphaFdd;
    const double CF = 1.0 + betaFddd*(5.0 + 4.0*betaFddd + 60.0*betaFddd*betaFddd) 
                    + 5.0*(1.0 + 3.0*betaFddd + 10.0*betaFddd*betaFddd)*alphaFddd 
                    + 2.0*(4.0 + 11.0*betaFddd)*alphaFddd*alphaFddd 
                    + 5.0*alphaFddd*alphaFddd*alphaFddd;

    const double yF00=0.0;
    const double yF10=(10.0*betaFdd*betaFdd*(8.0*betaFdd - 1.0) 
                     + (1.0 + 4.0*betaFdd + 81.0*betaFdd*betaFdd)*alphaFdd 
                     + 5.0*(1.0 + 8.0*betaFdd)*alphaFdd*alphaFdd 
                     + 9.0*alphaFdd*alphaFdd*alphaFdd)/BF;
    const double yF20=betaFd;
    const double yF01=(alphaFddd*(1.0 + alphaFddd)*(1.0 + 4.0*alphaFddd) 
                     + 2.0*alphaFddd*(7.0 + 3.0*alphaFddd)*betaFddd 
                     + 24.0*(1.0 - alphaFddd)*betaFddd*betaFddd 
                     - 80.0*betaFddd*betaFddd*betaFddd)/CF;
    const double yF11=0.0;
    const double yF21=alphaFd;
    const double yF02=(alphaFddd*alphaFddd*alphaFddd 
                     + (1.0 + 3.0*alphaFddd + 14.0*alphaFddd*alphaFddd)*betaFddd 
                     + 46.0*alphaFddd*betaFddd*betaFddd 
                     + 60.0*betaFddd*betaFddd*betaFddd)/CF;
    const double yF12=(alphaFdd*(1.0 + 5.0*alphaFdd + 9.0*alphaFdd*alphaFdd) 
                     + alphaFdd*(5.0 + 36.0*alphaFdd)*betaFdd 
                     + (55.0*alphaFdd - 1.0)*betaFdd*betaFdd 
                     + 10.0*betaFdd*betaFdd*betaFdd) / BF;
    const double yF22=0.0;
    const double yF03=0.0;
    const double yF13=betaFdd*(1.0 + 5.0*alphaFdd + 9.0*alphaFdd*alphaFdd + 5.0*(1.0 + 7.0*alphaFdd)*betaFdd + 50.0*betaFdd*betaFdd)/BF;
    const double yF23=alphaFd;
    const double yF04=0.0;
    const double yF14=0.0;
    const double yF24=betaFd;

    const double bF20=aF2d +  5.0*aF3d;
    const double bF21=aF1d - 10.0*aF3d;
    const double bF23=aF1d -  5.0*aF3d;
    const double bF24=aF2d + aF3d;
    const double bF25=aF3d;
    const double bF22= -(bF20 + bF21 + bF23 + bF24 + bF25);



    const int nd = n * n;

    /* initialize the matrix to zero. */
    for (int i = 0; i < nd; i++) {
        P[i] = 0.0;
    }

    for (int i = 3; i < n-3; i++) {
        P[IDX(i,i-2)] = betaF;
        P[IDX(i,i-1)] = alphaF;
        P[IDX(i,i)] = 1.0;
        P[IDX(i,i+1)] = alphaF;
        P[IDX(i,i+2)] = betaF;
    }

    P[IDX(0,0)] = 1.0;
    P[IDX(0,1)] = yF01;
    P[IDX(0,2)] = yF02;

    P[IDX(1,0)] = yF10;
    P[IDX(1,1)] = 1.0;
    P[IDX(1,2)] = yF12;
    P[IDX(1,3)] = yF13;

    P[IDX(2,0)] = yF20;
    P[IDX(2,1)] = yF21;
    P[IDX(2,2)] = 1.0;
    P[IDX(2,3)] = yF23;
    P[IDX(2,4)] = yF24;

    P[IDX(n-3,n-5)] = yF24;
    P[IDX(n-3,n-4)] = yF23;
    P[IDX(n-3,n-3)] = 1.0;
    P[IDX(n-3,n-2)] = yF21;
    P[IDX(n-3,n-1)] = yF20;

    P[IDX(n-2,n-4)] = yF13;
    P[IDX(n-2,n-3)] = yF12;
    P[IDX(n-2,n-2)] = 1.0;
    P[IDX(n-2,n-1)] = yF10;

    P[IDX(n-1,n-3)] = yF02;
    P[IDX(n-1,n-2)] = yF01;
    P[IDX(n-1,n-1)] = 1.0;

    /* initialize the matrix to zero. */
    for (int i = 0; i < nd; i++) {
        Q[i] = 0.0;
    }

    for (int i = 3; i < n-3; i++) {
        Q[IDX(i,i-3)] = aF3;
        Q[IDX(i,i-2)] = aF2;
        Q[IDX(i,i-1)] = aF1;
        Q[IDX(i,i)]   = aF0;
        Q[IDX(i,i+1)] = aF1;
        Q[IDX(i,i+2)] = aF2;
        Q[IDX(i,i+3)] = aF3;
    }

    Q[IDX(0,0)] = 0.0;
    Q[IDX(0,1)] = 0.0;
    Q[IDX(0,2)] = 0.0;
    Q[IDX(0,3)] = 0.0;

    Q[IDX(1,0)] = 0.0;
    Q[IDX(1,1)] = 0.0;
    Q[IDX(1,2)] = 0.0;
    Q[IDX(1,3)] = 0.0;
    Q[IDX(1,4)] = 0.0;

    Q[IDX(2,0)] = bF20;
    Q[IDX(2,1)] = bF21;
    Q[IDX(2,2)] = bF22;
    Q[IDX(2,3)] = bF23;
    Q[IDX(2,4)] = bF24;
    Q[IDX(2,5)] = bF25;
 
    Q[IDX(n-3,n-6)] = bF25;
    Q[IDX(n-3,n-5)] = bF24;
    Q[IDX(n-3,n-4)] = bF23;
    Q[IDX(n-3,n-3)] = bF22;
    Q[IDX(n-3,n-2)] = bF21;
    Q[IDX(n-3,n-1)] = bF20;

    Q[IDX(n-2,n-5)] = 0.0;
    Q[IDX(n-2,n-4)] = 0.0;
    Q[IDX(n-2,n-3)] = 0.0;
    Q[IDX(n-2,n-2)] = 0.0;
    Q[IDX(n-2,n-1)] = 0.0;

    Q[IDX(n-1,n-4)] = 0.0;
    Q[IDX(n-1,n-3)] = 0.0;
    Q[IDX(n-1,n-2)] = 0.0;
    Q[IDX(n-1,n-1)] = 0.0;

}
#endif // KIM_FILTER_H
// }}}