#include "transform.h"

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

    for (i = n_points; i > 0; i--) {
        fscanf (f, "%i,%lf", &dum, &data);
        GSL_SET_COMPLEX(&z, data, 0);
        gsl_vector_complex_set (dest, i-1, z);
    }
}

void printArray(FILE *f, gsl_vector_complex *source, int n_points) {
    int i;
    for (i = 0; i < n_points; i++)
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
    for (i = 0; i < n_points; i++) {
        a_val = GSL_REAL(gsl_vector_complex_get(a,i));
        b_val = GSL_REAL(gsl_vector_complex_get(b,i));
        sum += (a_val - b_val)*(a_val - b_val);
    }
    return sqrt(sum/n_points);
}

int main(int argc, char **argv) {
    int i, n_points;
    char o_name[] = "out";
    FILE *input, *out, **dtft_out, **dft_out;
    char buffer[100];
    double time, error;
    gsl_vector_complex *raw;
    gsl_vector_complex *trans;

    if (argc != 3) {
        printf ("Usage: %s LINES FILE\n", argv[0]);
        exit(1);
    }

    /* Get commandline args */
    n_points = atoi(argv[1]);

    /* Allocate memory */
    //puts("Allocating memory.");
    dft_out = malloc(4 * sizeof(FILE*));
    dtft_out = malloc(4 * sizeof(FILE*));
    raw = gsl_vector_complex_alloc(n_points);
    trans = gsl_vector_complex_alloc(n_points);
    if (dft_out == NULL || dtft_out == NULL) {
        printf ("Couldn't allocate memory.\n");
        exit(1);
    }

    /* Open files */
    //puts("Opening files.");
    input = my_fopen(argv[2], "r");
    out = my_fopen(o_name, "w");
    strcpy(buffer, o_name);
    for (i = 0; i < 4; i++) {
        sprintf(buffer, "%s.dft.%i", o_name, 2*i+1);
        dft_out[i] = my_fopen(buffer, "w");
        sprintf(buffer, "%s.dtft.%i", o_name, 2*i+1);
        dtft_out[i] = my_fopen(buffer, "w");
    }

    /* Read in the raw data */
    //puts("Reading in data.");
    getArrayFromFile(input, raw, n_points);
    printArray(out, raw, n_points);
    //gsl_vector_complex_fprintf(stdout, raw, "%e");

    /* Make the discrete tranforms */
    for (i = 0; i < 4; i++) {
        //printf ("Transforming with %i frequencies.\n", 2*i+1);
        if (i == 0)
            time = getDFTFromArray(trans, raw, n_points, 2*i+1, TRUE);
        else
            time = getDFTFromArray(trans, NULL, n_points, 2*i+1, TRUE);
        error = getRMSError(raw, trans, n_points);
        printf ("Got DFT with %i frequencies in %lf secs with error %e.\n", 2*i+1, time, error);
        printArray(dft_out[i], trans, n_points);
    }

    /* Make the discrete time series */
    for (i = 0; i < 4; i++) {
        time = getDFTFromArray(trans, NULL, n_points, 2*i+1, FALSE);
        error = getRMSError(raw, trans, n_points);
        printf ("Got DTFT with %i frequencies in %lf secs with error %e.\n", 2*i+1, time, error);
        printArray(dtft_out[i], trans, n_points);
    }
    getDFTFromArray(NULL, NULL, n_points, 2*i+1, 0);

    /* Close all files and free memory */
    fclose(input);
    fclose(out);
    for (i = 0; i < 4; i++) {
        fclose(dft_out[i]);
        fclose(dtft_out[i]);
    }

    exit(0);
}
