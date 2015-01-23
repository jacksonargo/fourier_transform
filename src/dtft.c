/* src/dtft.c */

#ifdef DONOTCOMPILE
#include "transform.h"

gsl_complex getSeriesCoeffs(gsl_vector_complex *v, gsl_maxtrix_complex *FOR, int n_points, int k) {
    int i;
    gsl_complex z, w, sum;

    for (i = 0; i < n_points; i++) {
        w = gsl_vector_complex_get(v, i);
        GSL_SET_COMPLEX(&z, 0.0, -2.0*pi*k*i/n_points);
        z = gsl_complex_exp(z);
        z = gsl_complex_mul(z,w);
        sum = gsl_complex_add(sum,z);
    }
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, , u, GSL_COMPLEX_ZERO, dest);
    return ;
}

double getDTFTFromArray(gsl_vector_complex *dest, gsl_vector_complex *source, gsl_matrix_complex *INV, int n_points, int n_frequencies) {
    static gsl_vector_complex *c;
    static int prev_n_freqs;
    gsl_complex z;
    int i;


    if (dest == NULL || source != NULL) {
        gsl_vector_complex_free(c);
        prev_n_freqs = 0;

        if (dest == NULL)
            return 0.0;

        c = gsl_vector_complex_calloc (n_points);
    }

    /* Add or remove from coeff list */

    if (prev_n_freqs < n_frequencies) {
        for (i = prev_n_freqs; i < n_frequencies; i++) {
            z = getSeriesCoeffs(source, n_points, i);
            gsl_vector_complex_set(c, i, z);
        }
    }

    else if (prev_n_freqs > n_frequencies) {
        for (i = n_frequencies; i < prev_n_freqs; i++) {
            gsl_vector_complex_set(c, i, GSL_COMPLEX_ZERO);
        }
    }

    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, INV, c, GSL_COMPLEX_ZERO, dest);
    prev_n_freqs = n_frequencies;
    return timer();
}
#endif
