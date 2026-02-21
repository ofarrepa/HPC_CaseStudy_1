#include "tsqr.h"
#include <stdlib.h>
#include <string.h>

extern void dgeqrf_(int* M, int* N, double* A, int* LDA, double* TAU,
                    double* WORK, int* LWORK, int* INFO);
extern void dorgqr_(int* M, int* N, int* K, double* A, int* LDA, double* TAU,
                    double* WORK, int* LWORK, int* INFO);
extern void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
                   double* ALPHA, double* A, int* LDA, double* B, int* LDB,
                   double* BETA, double* C, int* LDC);

static void extract_R(const double* Aqr, int lda, int n, double* R) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            R[i + j*n] = (i <= j) ? Aqr[i + j*lda] : 0.0;
        }
    }
}

static void copy_block_rows(const double* A, int m, int n, int r0, int mb, double* B) {
    for (int j = 0; j < n; j++) {
        memcpy(B + (size_t)j*mb,
               A + (size_t)j*m + r0,
               (size_t)mb * sizeof(double));
    }
}


static void place_block_rows(double* Q, int m, int n, int r0, int mb, const double* C) {
    for (int j = 0; j < n; j++) {
        memcpy(Q + (size_t)j*m + r0,
               C + (size_t)j*mb,
               (size_t)mb * sizeof(double));
    }
}

static int local_qr_explicitQ(const double* Ablk, int mb, int n,
                             double* Qblk, double* Rblk) {
    memcpy(Qblk, Ablk, (size_t)mb * (size_t)n * sizeof(double));

    int info = 0;
    int M = mb, N = n, LDA = mb;

    double* tau = (double*)malloc((size_t)n * sizeof(double));
    if (!tau) return -1;

    //dgeqrf quer//
    int lwork = -1;
    double wkopt = 0.0;
    dgeqrf_(&M, &N, Qblk, &LDA, tau, &wkopt, &lwork, &info);
    if (info != 0) { free(tau); return -2; }

    lwork = (int)wkopt;
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) { free(tau); return -1; }

    //the factor//
    dgeqrf_(&M, &N, Qblk, &LDA, tau, work, &lwork, &info);
    if (info != 0) { free(work); free(tau); return -2; }

    extract_R(Qblk, LDA, n, Rblk);

    //dorgqr query again//
    lwork = -1;
    dorgqr_(&M, &N, &N, Qblk, &LDA, tau, &wkopt, &lwork, &info);
    if (info != 0) { free(work); free(tau); return -3; }

    lwork = (int)wkopt;
    double* work2 = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work2) { free(work); free(tau); return -1; }

    //forming th explicit Q//
    dorgqr_(&M, &N, &N, Qblk, &LDA, tau, work2, &lwork, &info);
    if (info != 0) { free(work2); free(work); free(tau); return -3; }

    free(work2);
    free(work);
    free(tau);
    return 0;
}

int TSQR(const double* A, int m, int n, double* Q, double* R) {
    if (m < n) return -10;

    //Splitting into 4row blocks//
    int c0 = 0;
    int c1 = m / 4;
    int c2 = m / 2;
    int c3 = (3 * m) / 4;
    int c4 = m;

    int m1 = c1 - c0;
    int m2 = c2 - c1;
    int m3 = c3 - c2;
    int m4 = c4 - c3;

    //allocating the AQR blocks// 
    double *A1 = (double*)malloc((size_t)m1*(size_t)n*sizeof(double));
    double *A2 = (double*)malloc((size_t)m2*(size_t)n*sizeof(double));
    double *A3 = (double*)malloc((size_t)m3*(size_t)n*sizeof(double));
    double *A4 = (double*)malloc((size_t)m4*(size_t)n*sizeof(double));

    double *Q1 = (double*)malloc((size_t)m1*(size_t)n*sizeof(double));
    double *Q2 = (double*)malloc((size_t)m2*(size_t)n*sizeof(double));
    double *Q3 = (double*)malloc((size_t)m3*(size_t)n*sizeof(double));
    double *Q4 = (double*)malloc((size_t)m4*(size_t)n*sizeof(double));

    double *R1 = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    double *R2 = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    double *R3 = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    double *R4 = (double*)malloc((size_t)n*(size_t)n*sizeof(double));

    if (!A1||!A2||!A3||!A4||!Q1||!Q2||!Q3||!Q4||!R1||!R2||!R3||!R4) return -1;

    copy_block_rows(A, m, n, c0, m1, A1);
    copy_block_rows(A, m, n, c1, m2, A2);
    copy_block_rows(A, m, n, c2, m3, A3);
    copy_block_rows(A, m, n, c3, m4, A4);

    int rc = 0;
    rc = local_qr_explicitQ(A1, m1, n, Q1, R1); if (rc) return rc;
    rc = local_qr_explicitQ(A2, m2, n, Q2, R2); if (rc) return rc;
    rc = local_qr_explicitQ(A3, m3, n, Q3, R3); if (rc) return rc;
    rc = local_qr_explicitQ(A4, m4, n, Q4, R4); if (rc) return rc;

    //Stacking the R's into Rstack (4n x n)//
    int M = 4*n;
    int N = n;
    int LDA = M;

    double* Rstack = (double*)calloc((size_t)M*(size_t)N, sizeof(double));
    if (!Rstack) return -1;

    //column-by-column stacking into Rstack which is column major// 
    const double* Rblk[4] = {R1, R2, R3, R4};
    for (int bi = 0; bi < 4; ++bi) {
        int row0 = bi * n;
        for (int j = 0; j < n; ++j) {
            memcpy(Rstack + row0 + (size_t)j*M,
                   Rblk[bi] + (size_t)j*n,
                   (size_t)n * sizeof(double));
        }
    }

    //QR of Rstack which then stored in QR, then forms explicit Q_R and final R//
    double* QR = (double*)malloc((size_t)M*(size_t)N*sizeof(double));
    if (!QR) return -1;
    memcpy(QR, Rstack, (size_t)M*(size_t)N*sizeof(double));

    double* tau = (double*)malloc((size_t)n * sizeof(double));
    if (!tau) return -1;

    int info = 0;
    
    int lwork = -1;
    double wkopt = 0.0;
    dgeqrf_(&M, &N, QR, &LDA, tau, &wkopt, &lwork, &info);
    if (info != 0) return -2;

    lwork = (int)wkopt;
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -1;

    dgeqrf_(&M, &N, QR, &LDA, tau, work, &lwork, &info);
    if (info != 0) return -2;

    extract_R(QR, LDA, n, R);

    //forming explicit Q_R in QR
    lwork = -1;
    dorgqr_(&M, &N, &N, QR, &LDA, tau, &wkopt, &lwork, &info);
    if (info != 0) return -3;

    lwork = (int)wkopt;
    double* work2 = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work2) return -1;

    dorgqr_(&M, &N, &N, QR, &LDA, tau, work2, &lwork, &info);
    if (info != 0) return -3;

    //assembling the final Q//
    char transN = 'N';
    double alpha = 1.0, beta = 0.0;
    int K = n;
    int ldb = M; 

    double *C1 = (double*)malloc((size_t)m1*(size_t)n*sizeof(double));
    double *C2 = (double*)malloc((size_t)m2*(size_t)n*sizeof(double));
    double *C3 = (double*)malloc((size_t)m3*(size_t)n*sizeof(double));
    double *C4 = (double*)malloc((size_t)m4*(size_t)n*sizeof(double));
    if (!C1||!C2||!C3||!C4) return -1;

    int lda, ldc;

    lda = m1; ldc = m1;
    dgemm_(&transN, &transN, &m1, &n, &K, &alpha, Q1, &lda, QR + 0*n, &ldb, &beta, C1, &ldc);

    lda = m2; ldc = m2;
    dgemm_(&transN, &transN, &m2, &n, &K, &alpha, Q2, &lda, QR + 1*n, &ldb, &beta, C2, &ldc);

    lda = m3; ldc = m3;
    dgemm_(&transN, &transN, &m3, &n, &K, &alpha, Q3, &lda, QR + 2*n, &ldb, &beta, C3, &ldc);

    lda = m4; ldc = m4;
    dgemm_(&transN, &transN, &m4, &n, &K, &alpha, Q4, &lda, QR + 3*n, &ldb, &beta, C4, &ldc);

    place_block_rows(Q, m, n, c0, m1, C1);
    place_block_rows(Q, m, n, c1, m2, C2);
    place_block_rows(Q, m, n, c2, m3, C3);
    place_block_rows(Q, m, n, c3, m4, C4);

    //freeing memory// 
    free(C4); free(C3); free(C2); free(C1);
    free(work2); free(work); free(tau);
    free(QR); free(Rstack);
    free(R4); free(R3); free(R2); free(R1);
    free(Q4); free(Q3); free(Q2); free(Q1);
    free(A4); free(A3); free(A2); free(A1);

    return 0;
}
