#include "tsqr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <time.h> 

static double wall_seconds(void) { 
	 return(double)clock() / (double)CLOCKS_PER_SEC; 
} 



static double randu(void) {
    return (rand() + 1.0) / (RAND_MAX + 2.0);
}

static double randn(void) {
    double u1 = randu(), u2 = randu();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void fill_randn(double* A, int m, int n) {
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            A[i + j*m] = randn();
}

static double best_time_tsqr(int m, int n, int reps) {
    double* A = (double*)malloc((size_t)m*(size_t)n*sizeof(double));
    double* Q = (double*)malloc((size_t)m*(size_t)n*sizeof(double));
    double* R = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    if (!A || !Q || !R) { fprintf(stderr, "alloc failed\n"); exit(1); }

    fill_randn(A, m, n);

    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        double t0 = wall_seconds();
        TSQR(A, m, n, Q, R);
        double t1 = wall_seconds();
        double dt = t1 - t0;
        if (dt < best) best = dt;
    }

    free(R); free(Q); free(A);
    return best;
}

int main(void) {
    srand(0);

    FILE* f = fopen("scaling.csv", "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "sweep,m,n,seconds\n");

    int n_fixed = 64;
    int ms[] = {2000, 5000, 10000, 20000, 40000, 80000};
    int nm = (int)(sizeof(ms)/sizeof(ms[0]));

    for (int i = 0; i < nm; ++i) {
        int m = ms[i];
        double t = best_time_tsqr(m, n_fixed, 3);
        fprintf(f, "m_sweep,%d,%d,%.9f\n", m, n_fixed, t);
        fflush(f);
        printf("m_sweep m=%d n=%d time=%.6f s\n", m, n_fixed, t);
    }

    int m_fixed = 40000;
    int ns[] = {16, 32, 64, 96, 128, 192, 256};
    int nn = (int)(sizeof(ns)/sizeof(ns[0]));

    for (int i = 0; i < nn; ++i) {
        int n = ns[i];
        double t = best_time_tsqr(m_fixed, n, 3);
        fprintf(f, "n_sweep,%d,%d,%.9f\n", m_fixed, n, t);
        fflush(f);
        printf("n_sweep m=%d n=%d time=%.6f s\n", m_fixed, n, t);
    }

    fclose(f);
    printf("Wrote scaling.csv\n");
    return 0;
}
