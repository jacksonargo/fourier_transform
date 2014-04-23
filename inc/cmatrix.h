/* inc/cmatrix.h */

#ifndef __CMATRIX__
#define __CMATRIX__

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <float.h>

typedef struct __cRow cRow;
struct __cRow {
    double complex *col;
};

typedef struct __cmatrix cmatrix;
struct __cmatrix {
    cRow *row;
    int n;
    int m;
};

cmatrix cinitmatrix (int n, int m);
cmatrix cmatrixmultiply (cmatrix A, cmatrix B, cmatrix *S);
void cprintmatrix (cmatrix A);
cmatrix cscalematrix (cmatrix A, double c);
cmatrix conjugatematrix (cmatrix A);
cmatrix cduplicatematrix (cmatrix A);

#endif
