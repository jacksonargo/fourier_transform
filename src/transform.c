#include "transform.h"

cmatrix makeTransMatrix (int n) {
    int i, j;
    cmatrix F = cinitmatrix(n,n);
    for (i = 0; i < F.n; i++)
        for (j = 0; j < F.m; j++)
            F.row[i].col[j] = cexp (-2*pi*I/n*i*j);
    return F;
}

cmatrix makeRandomVector(int n) {
    cmatrix A = cinitmatrix(n,1);
    int i;
    for (i = 0; i < A.n; i++)
        A.row[i].col[0] = RAND_DBL*sin((double)i*6*pi/A.n);
    return A;
}

void printRealToStream(FILE *stream, cmatrix A) {
    int i, j;
    for (i = 0; i < A.n; i++) {
        for (j = 0; j < A.m; j++) {
            fprintf (stream, "%.*e ", DBL_DIG, creal(A.row[i].col[j]));
        }
        fprintf (stream, "\n");
    }
}

int cmpDouble(const void *pa, const void *pb) {
    double a = *(double*)pa;
    double b = *(double*)pb;

    if (a > b)
        return 1;
    else if (a < b)
        return -1;
    else
        return 0;
}

int cmpVectorElems(const void *pa, const void *pb) {
    double a = cabs(((cRow*)pa)->col[0]);
    double b = cabs(((cRow*)pb)->col[0]);

    return cmpDouble(&a, &b);
}

cmatrix sortVector (cmatrix A) {
    int i;
    cmatrix S;

    S.n = A.n;
    S.m = 1;
    S.row = malloc(S.n*sizeof(cRow*));
    for (i = 0; i < S.n; i++) {
        S.row[i] = A.row[i];
    }

    qsort(S.row, A.n, sizeof(cRow), cmpVectorElems);
    return S;
}

void freeRowsCols (void *pA) {
    cmatrix A = *(cmatrix*)pA;
    int i;
    for (i = 0; i < A.n; i++)
        free(A.row[i].col);
    free(A.row);
}

double getError(cmatrix A, cmatrix B) {
    int i;
    double error = 0;
    for (i = 0; i < A.n; i++)
        error += (A.row[i].col[0]-B.row[i].col[0])*(A.row[i].col[0]-B.row[i].col[0]);
    return sqrt(fabs(error)/A.n);
}

cmatrix removeFrequencies (cmatrix s, double allowable_error) {
    cmatrix F = makeTransMatrix(s.n); // Make the DSF matrix and inverse
    cmatrix G = cscalematrix(conjugatematrix(cduplicatematrix(F)), 1.0/s.n);
    cmatrix t = cmatrixmultiply(F,s,NULL); // Transform s into frequency domain
    cmatrix w = cduplicatematrix(t);  // Copy it to w
    cmatrix u = sortVector(w);        // Sort w
    cmatrix v = cinitmatrix(s.n,1);
    double error;
    int i, j;
    FILE *p = fopen("error", "w");

    // Delete all but the first element
    j = 1;
    for(;;) {
        for (i = 0; i < s.n - j; i++)
            u.row[i].col[0] = 0;
        v = cmatrixmultiply(G,w,&v); // Convert back to time domain
        error = getError(s,v); // Check the error
        fprintf(p, "%i %.*e\n", j, DBL_DIG, error);
        //if (error > allowable_error) { // If the error is too big
        if (j < s.n) { // I want error at each spot
            j++; // Increase the number we keep
            //freeRowsCols(&v); // Free v
            for (i = 0; i < s.n; i++) // Reset w back to t
                w.row[i].col[0] = t.row[i].col[0];
        }
        else // We're done
            break;
    }

    freeRowsCols(&F);
    freeRowsCols(&G);
    freeRowsCols(&t);
    printf ("Error with %i frequencies: %.*e.\n", j, DBL_DIG, error);
    fclose(p);
    return v;
}

int main(int argc, char **argv) {
    //srand(time(NULL));
    double r = (argc > 1) ? atof(argv[1]) : 1;
    int n = 1000;
    FILE *a = fopen ("original", "w");
    FILE *b = fopen ("transformed", "w");
    cmatrix s = makeRandomVector(n);
    cmatrix t = removeFrequencies(s, r);

    printRealToStream(a, s);
    printRealToStream(b, t);

    fclose(a);
    fclose(b);
    exit(0);
}

