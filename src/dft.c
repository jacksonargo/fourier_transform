/* src/dft.c */

#include "transform.h"

void makeFourierMatrix(gsl_matrix_complex *F, int n_points) {
    int i,j;
    gsl_complex z;
    for (i = 0; i < n_points; i++) {
        for (j = 0; j < n_points; j++) {
            GSL_SET_COMPLEX(&z, 0.0, -2.0*pi*i*j/n_points);
            z = gsl_complex_exp(z);
            //printf("%e %e\n", GSL_REAL(z), GSL_IMAG(z));
            gsl_matrix_complex_set(F, i, j, z);
        }
    }
}

void makeInverseMatrix(gsl_matrix_complex *G, gsl_matrix_complex *F, int n_points) {
    int i, j;
    gsl_complex z;
    for (i = 0; i < n_points; i++) {
        for (j = 0; j < n_points; j++) {
            z = gsl_matrix_complex_get(F, i, j);       // Copy
            z = gsl_complex_conjugate(z);              // Conjugate
            z = gsl_complex_div_real(z, (double)n_points); // Scale
            gsl_matrix_complex_set(G, i, j, z);
        }
    }
}

double getDFTFromArray(gsl_vector_complex *dest, gsl_vector_complex *source, int n_points, int n_frequencies) {
    static int first_run = TRUE;
    static gsl_matrix_complex *F, *G; 
    static gsl_vector_complex *v, *u;
    static gsl_vector *m;
    static gsl_permutation *p;
    double mag;
    int i, j;

    if (first_run) {
        //puts("First call of getDFTFromArray so allocating memory");
        first_run = FALSE;
        F = gsl_matrix_complex_alloc (n_points, n_points);
        G = gsl_matrix_complex_alloc (n_points, n_points);
        v = gsl_vector_complex_alloc (n_points);
        u = gsl_vector_complex_alloc (n_points);
        m = gsl_vector_alloc (n_points);
        p = gsl_permutation_alloc(n_points);

        timer(); // Start the timer.
        //puts("Making Fourier matrix.");
        makeFourierMatrix(F, n_points);
        //gsl_matrix_complex_fprintf(stdout, F, "%e");
        //puts("Making Inverse Fourier matrix.");
        makeInverseMatrix(G, F, n_points);
        //gsl_matrix_complex_fprintf(stdout, G, "%e");



        /* Transform source */
        //puts("Transforming the source");
        gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, F, source, GSL_COMPLEX_ZERO, v);
        //gsl_vector_complex_fprintf(stdout, v, "%e");

        /* Store the magnitudes in m */
        //puts("Getting and storing magnitudes.");
        for (i = 0; i < n_points; i++) {
            mag = gsl_complex_abs(gsl_vector_complex_get(v, i));
            gsl_vector_set (m, i, mag);
        }
        //gsl_vector_fprintf(stdout, m, "%e");
        /* Sort m */
        //puts("Sorting magnitudes.");
        gsl_sort_vector_index(p, m);
    }
    else
        timer(); // Start the timer.

    /* We need to duplicate v */
    //puts("Making duplicate of transformed vector.");
    gsl_vector_complex_memcpy(u, v);
    //puts("");gsl_vector_complex_fprintf(stdout, u, "%e");

    /* Remove frequencies */
    //printf ("Removing frequencies.\n");
    for (i = 0; i < n_points - n_frequencies; i++) {
        j = gsl_permutation_get (p, i);
        gsl_vector_complex_set (u, j, GSL_COMPLEX_ZERO);
    }
    //puts("");gsl_vector_complex_fprintf(stdout, u, "%e");

    /* Tranform back to time */
    //puts("Transforming back");
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, G, u, GSL_COMPLEX_ZERO, dest);
    //puts("");gsl_vector_complex_fprintf(stdout, dest, "%e");

    return timer();
}
