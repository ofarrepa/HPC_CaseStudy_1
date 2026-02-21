#include "tsqr.h" 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static double randu(void) {
    return (rand()+1.0)/(RAND_MAX+2.0);
}

static double randn(void) {
    double u1 = randu(), u2 = randu();
    return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

static double frob(const double* A, int m, int n) {
    double s=0;
    for(int j=0;j<n;j++)
        for(int i=0;i<m;i++) s+=A[i+j*m]*A[i+j*m];
    return sqrt(s);
}

int main(void) {
    srand(0);
    int m=2000, n=40;

    double *A=malloc(m*n*sizeof(double));
    double *Q=malloc(m*n*sizeof(double));
    double *R=malloc(n*n*sizeof(double));
    double *QR=malloc(m*n*sizeof(double));

    for(int j=0;j<n;j++)
        for(int i=0;i<m;i++)
            A[i+j*m]=randn();

    TSQR(A,m,n,Q,R);

    memset(QR,0,m*n*sizeof(double));
    for(int j=0;j<n;j++)
        for(int k=0;k<n;k++)
            for(int i=0;i<m;i++)
                QR[i+j*m]+=Q[i+k*m]*R[k+j*n];

    for(int i=0;i<m*n;i++) QR[i]=A[i]-QR[i];

    printf("Relative error = %.3e\n",frob(QR,m,n)/frob(A,m,n));
    return 0;
}
