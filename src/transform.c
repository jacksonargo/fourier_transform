#include "transform.h"

void fakeFillData(gsl_vector_complex *dest) {
    int n_head = 1024;
    int n_tail = 1024;
    int i;
    gsl_complex z,w;

    for (i = 0; i < n_head; i++) {
        z = gsl_vector_complex_get(dest, 1024);
        z = gsl_complex_mul_real(z,2);
        w = gsl_vector_complex_get(dest,1024+i);
        z = gsl_complex_sub(z,w);
        gsl_vector_complex_set(dest,1024-1-i,z);
    }
    for (i = 0; i < n_tail; i++) {
        z = gsl_vector_complex_get(dest, 3072);
        z = gsl_complex_mul_real(z,2);
        w = gsl_vector_complex_get(dest,3072-1-i);
        z = gsl_complex_sub(z,w);
        gsl_vector_complex_set(dest,i+3073,z);
    }
}
double timer () {
    static clock_t start = 0;
    static clock_t stop;
    double change;

    stop = clock();
    change = (double)(stop - start)/CLOCKS_PER_SEC;
    start = clock();
    return change;
}

void getArrayFromFile(FILE *f, gsl_vector_complex *dest, int n_points) {
    int i, dum;
    double data;
    gsl_complex z;

    for (i = 1024*3+1; i > 1024; i--) {
        fscanf (f, "%i,%lf", &dum, &data);
        GSL_SET_COMPLEX(&z, data, 0);
        gsl_vector_complex_set (dest, i-1, z);
    }

    fakeFillData(dest);
}

void getRombArray(FILE *f, gsl_vector_complex *dest, int n_points) {
    int i;
    double data;
    gsl_complex z;
    for (i = n_points; i > 0; i--) {
        fscanf (f, "%lf", &data);
        GSL_SET_COMPLEX(&z, data, 0);
        gsl_vector_complex_set (dest, i-1, z);
    }
}


void printArray(FILE *f, gsl_vector_complex *source, int n_points) {
    int i;
    //for (i = 0; i < n_points; i++)
    for (i = 1024; i < 3*1024+1; i++)
        fprintf (f, "%.*e\n", DBL_DIG, GSL_REAL(gsl_vector_complex_get(source, i)));
}

FILE *my_fopen(char *path, char *mode) {
    FILE *f = fopen(path, mode);
    if (f == NULL) {
        printf ("Couldn't access %s for reading.\n", path);
        exit(1);
    }
    return f;
}

double getRMSError(gsl_vector_complex *a, gsl_vector_complex *b, int n_points) {
    double a_val, b_val, sum = 0;
    int i;
    for (i = 1024; i < 1024*3+1; i++) {
        a_val = GSL_REAL(gsl_vector_complex_get(a,i));
        b_val = GSL_REAL(gsl_vector_complex_get(b,i));
        sum += (a_val - b_val)*(a_val - b_val);
    }
    return sqrt(sum/n_points);
}


int main(int argc, char **argv) {
    int i, n_points;
    char o_name[] = "out";
    FILE *input, *out, **dtft_out, **dft_out, *romb_in;
    char buffer[100];
    double time, error;
    gsl_vector_complex *raw, *romb;
    gsl_vector_complex *trans;

    /* Get commandline args */
    n_points = 4097;

    /* Allocate memory */
    //puts("Allocating memory.");
    dft_out = malloc(12 * sizeof(FILE*));
    dtft_out = malloc(12 * sizeof(FILE*));
    raw = gsl_vector_complex_alloc(n_points);
    romb = gsl_vector_complex_alloc(n_points);
    trans = gsl_vector_complex_alloc(n_points);
    if (dft_out == NULL || dtft_out == NULL) {
        printf ("Couldn't allocate memory.\n");
        exit(1);
    }

    /* Open files */
    //puts("Opening files.");
    romb_in = my_fopen("romb_in", "r");
    input = my_fopen("data.csv", "r");
    out = my_fopen(o_name, "w");
    strcpy(buffer, o_name);
    for (i = 0; i < 12; i++) {
        sprintf(buffer, "%s.dft.%i", o_name, 2*i+1);
        dft_out[i] = my_fopen(buffer, "w");
        sprintf(buffer, "%s.dtft.%i", o_name, 2*i+1);
        dtft_out[i] = my_fopen(buffer, "w");
    }

    /* Read in the raw data */
    //puts("Reading in data.");
    getArrayFromFile(input, raw, n_points);
    getRombArray(romb_in, romb, n_points);
    gsl_vector_complex_mul(romb, raw);
    printArray(out, raw, n_points);
    //printArray(out,romb,n_points);
    //puts("Done.");exit(0);
    //gsl_vector_complex_fprintf(stdout, raw, "%e");

    /* Make the discrete tranforms */
#ifdef DONOTCOMPILE
    for (i = 0; i < 4097; i++) {
        gsl_complex z;
        z = gsl_vector_complex_get(raw,i);
        z = gsl_complex_mul_real(z,1.0/3.0);
        if (i == 0 || i == 4096)
             z = gsl_complex_mul_real(z,1);
        else if (i%2 == 1)
             z = gsl_complex_mul_real(z,4);
        else if (i%2 == 0)
             z = gsl_complex_mul_real(z,2);
        gsl_vector_complex_set(raw,i,z);
    }
#endif
    romb = NULL;
    for (i = 0; i<12; i++) {
        //printf ("Transforming with %i frequencies.\n", 2*i+1);
        if (i == 0)
            time = getDFTFromArray(trans, raw, romb, n_points, 2*i+1, TRUE);
        else
            time = getDFTFromArray(trans, NULL, romb, n_points, 2*i+1, TRUE);
        error = getRMSError(raw, trans, n_points);
        printf ("Got DFT with %i frequencies in %lf secs with error %e.\n", 2*i+1, time, error);
        if (i < 12)
            printArray(dft_out[i], trans, n_points);
    }

    /* Make the discrete time series */
    for (i = 0; i < 12; i++) {
        time = getDFTFromArray(trans, NULL, romb, n_points, 2*i+1, FALSE);
        error = getRMSError(raw, trans, n_points);
        printf ("Got DTFT with %i frequencies in %lf secs with error %e.\n", 2*i+1, time, error);
        if (i < 12)
            printArray(dtft_out[i], trans, n_points);
    }
    getDFTFromArray(NULL, NULL, NULL, n_points, 2*i+1, 0);

    /* Close all files and free memory */
    fclose(input);
    fclose(out);
    for (i = 0; i < 12; i++) {
        fclose(dft_out[i]);
        fclose(dtft_out[i]);
    }

    exit(0);
}
