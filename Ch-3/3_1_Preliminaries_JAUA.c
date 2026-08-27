#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// De chat GPT sin depender del header nrutil.h
// Equivalente a interpolación de Polinomios de Laragnge??


/**
 * @brief Polynomial interpolation using Neville's algorithm.
 * 
 * Given arrays xa[0..n-1] and ya[0..n-1], and a value x, this function returns
 * an interpolated value y and an error estimate dy.
 * 
 * @param xa  Array of x data points (size n).
 * @param ya  Array of y data points (size n).
 * @param n   Number of data points.
 * @param x   Point at which interpolation is requested.
 * @param y   Pointer to store the output interpolated value.
 * @param dy  Pointer to store the output error estimate.
 * 
 * const make the dedicated variable read-only. It's value cannot be changed
 */

 
void polint(const double xa[], const double ya[], int n, double x, double *y, double *dy) {
    int i, m, ns = 0;       // index variables
                            // i -> find interval where the desired value x is located
                            // ns -> Index of the lower side of the desired interval. OptimiZes the calculatio to aviud not needed lower values
                            // m -> run over the table to calculate the desired polynomials.
                            //      i -> iruns over aech value of m
                            // 
    double den, dif, dift, ho, hp, w;

    /* Allocate working arrays for differences */
    double *c = (double *)malloc((size_t)n * sizeof(double));
    double *d = (double *)malloc((size_t)n * sizeof(double));

    if (c == NULL || d == NULL) {
        fprintf(stderr, "Allocation failure in polint\n");
        free(c);
        free(d);
        return;
    }

    /* Find the index of the closest table entry */
    dif = fabs(x - xa[0]);
    for (i = 0; i < n; i++) {
        if ((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        /* Initialize the C and D working arrays 
            Coefficints relain the value ya evaluated at xa
        */
        c[i] = ya[i];
        d[i] = ya[i];
    }

    /* Initial approximation to y */
    *y = ya[ns--];

    /* For each column of the tableau, loop over current c's and d's */
    for (m = 1; m < n; m++) {
        for (i = 0; i < n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];

            den = ho - hp;
            if (den == 0.0) { // Avoid division by zero
                fprintf(stderr, "Error in routine polint: two xa elements are equal!\n");
                free(c);
                free(d);
                return;
            }
            den = w / den;

            /* Update C and D */
            d[i] = hp * den;
            c[i] = ho * den;
        }

        /* 
         * Decide which correction (c or d) to add to y based on traversing
         * the straightest path through the tableau.
         */
        *dy = (2 * (ns + 1) < (n - m)) ? c[ns + 1] : d[ns--];
        *y += *dy;
    }

    free(c);
    free(d);
}



int main() {
    /* Sample data points: y = x^2 */
    double xa[] = {1.0, 2.0, 3.0, 4.0};
    double ya[] = {1.0, 4.0, 9.0, 16.0};
    int n = sizeof(xa) / sizeof(xa[0]);

    double x_test = 2.5;
    double y_result, dy_result;

    polint(xa, ya, n, x_test, &y_result, &dy_result);

    printf("Interpolating at x = %.2f:\n", x_test);
    printf("  y  = %.4f (Expected: %.4f)\n", y_result, x_test * x_test);
    printf("  dy = %.4e (Error estimate)\n", dy_result);

    return 0;
}