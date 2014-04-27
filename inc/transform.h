/* inc/transform.h */

#ifndef __TRANSFORM__
#define __TRANSFORM__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>

#define pi 3.1415926535897932384626433
#define RAND_DBL ((double)rand()/(RAND_MAX+1.0))

#define TRUE 1
#define FALSE 0

int main(int argc, char **argv);
double timer();

double getDFTFromArray(gsl_vector_complex *dest, gsl_vector_complex *source, int n_points, int n_frequencies, int maxed);
#endif
