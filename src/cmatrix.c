/* src/cmatrix.c */

#include "cmatrix.h"

cmatrix cinitmatrix (int n, int m) {
    int i;
    cmatrix F;
    F.n = n;
    F.m = m;
    F.row = calloc(F.n,sizeof(cRow*));
    for (i = 0; i < F.n; i++)
        F.row[i].col = calloc(F.m,sizeof(double complex));
    return F;
}

cmatrix cmatrixmultiply (cmatrix A, cmatrix B, cmatrix *S) {
    int i, j ,k;
    cmatrix C;
    if (S == NULL)
        C = cinitmatrix(A.n, B.m);
    else if ((*S).n != A.n || (*S).m != B.m)
        printf ("Matrix is not big enough to store product");
    else
        C = *S;

    if (A.m != B.n) {
        printf ("Incompatible sizes\n");
        exit(1);
    }

    for (i = 0; i < C.n; i++) {
        for (j = 0; j < C.m; j++) {
            C.row[i].col[j] = 0;
            for (k = 0; k < A.m; k++) {
                C.row[i].col[j] += A.row[i].col[k] * B.row[k].col[j];
            }
        }
    }
    return C;
}

void cprintmatrix (cmatrix A) {
    int i, j;
    for (i = 0; i < A.n; i++) {
        for (j = 0; j < A.m; j++) {
            printf ("%.*e%+.*ei ", DBL_DIG, creal(A.row[i].col[j]), DBL_DIG, cimag(A.row[i].col[j]));
        }
        printf ("\n");
    }
}

cmatrix cscalematrix (cmatrix A, double c) {
    int i, j;
    for (i = 0; i < A.n; i++)
        for (j = 0; j < A.m; j++)
            A.row[i].col[j] *= c;
    return A;
}

cmatrix conjugatematrix (cmatrix A) {
    int i, j;
    for (i = 0; i < A.n; i++)
        for (j = 0; j < A.m; j++)
            A.row[i].col[j] = conj(A.row[i].col[j]);
    return A;
}

cmatrix cduplicatematrix (cmatrix A) {
    int i, j;
    cmatrix C = cinitmatrix(A.n, A.m);
    for (i = 0; i < A.n; i++)
        for (j = 0; j < A.m; j++)
            C.row[i].col[j] = A.row[i].col[j];
    return C;
}

