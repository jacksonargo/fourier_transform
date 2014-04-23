/* inc/transform.h */

#ifndef __TRANSFORM__
#define __TRANSFORM__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>
#include "cmatrix.h"

#define pi 3.1415926535897932384626433
#define RAND_DBL ((double)rand()/(RAND_MAX+1.0))

cmatrix makeTransMatrix (int n);
cmatrix makeRandomVector(int n);
void printRealToStream(FILE *stream, cmatrix A);
int cmpDouble(const void *pa, const void *pb);
int cmpVectorElems(const void *pa, const void *pb);
cmatrix sortVector (cmatrix A);
void freeRowsCols (void *pA);
double getError(cmatrix A, cmatrix B);
cmatrix removeFrequencies (cmatrix s, double allowable_error);
int main(int argc, char **argv);

#endif
