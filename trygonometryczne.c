#include "makespl.h"
#include "gaus/matrix.h"
#include <math.h>
#include <stdlib.h>
#include "gaus/piv_ge_solver.h"

#ifndef APPROX_BASE_SIZE
#define APPROX_BASE_SIZE 4;
#endif

double diff_sin(double x, int n) {
    return sin(n * x);
}

double diff_cos(double x, int n) {
    return cos(n * x);
}

double *get_factors(matrix_t *);

void make_spl(points_t *pts, spline_t *spl) {
    int i, j, k;
    int parts = APPROX_BASE_SIZE;
    int x_size = pts->n;

    if (parts < x_size)
        fprintf(stderr, "\tIlość punktów jest większa niż funkcji bazowych - wynik znacząco traci na dokładności");

    matrix_t *equations = make_matrix(parts, parts + 1);

    /*zmienna k odpowiada za rząd - czyli dzięki k zmienia się: POCHODNA.*/
    for (k = 0; k < parts; k++) {
        double (*derivative)(double, int);

        /*n takie, że derivative=sin/cos(n*x)*/
        int derN = k / 2 + 1;

        /*jeśli k jest nieparzyste, to jesteśmy na sinusie (póki co - później dojdzie jeszcze liczba rzeczywista)*/
        if (k % 2 == 0) derivative = &diff_sin;
        else derivative = &diff_cos;

        /*zmienna j odpowiada za kolumnę - czyli zmienia: FUNKCJĘ.*/
        for (j = 0; j < parts + 1; j++) {
            /*suma w kolumnie*/
            double sum = 0;
            double (*func)(double);

            /*n takie, że func=sin/cos(n*x)*/
            double funcN = j / 2 + 1;

            if (j % 2 == 0) func = &sin;
            else func = &cos;

            if (j < parts) {
                //dodajemy pochodną razy wartość funkcji bazowej dla danego x
                for (i = 0; i < x_size; i++) {
                    double x = pts->x[i];

                    sum += func(x * funcN)
                           * derivative(x, derN);
                }
            }
            else {
                //dodajemy pochodną dla x razy wartość y które znamy
                printf("Pochodna: %s(%dx)\n", k % 2 ? "cos" : "sin", derN);
                for (i = 0; i < x_size; i++) {
                    double x = pts->x[i];
                    double y = pts->y[i];

                    sum += y * derivative(x, derN);
                }
            }

            put_entry_matrix(equations, k, j, sum);
        }
    }

    write_matrix(equations, stdout);

    piv_ge_solver(equations);

    write_matrix(equations, stdout);
    double *factors = get_factors(equations);
}

double *get_factors(matrix_t *mx) {
    int i;
    double *results = malloc(mx->rn * sizeof(double));
    for (i = 0; i < mx->rn; i++) {
        results[i] = get_entry_matrix(mx, i, mx->cn - 1);
    }
    return results;
}


//debugging
int main(int argc, char *argv[]) {
    points_t pts;
    spline_t spl;
    pts.n = 0;
    spl.n = 0;
    FILE *in = fopen("/Users/kubukoz/ClionProjects/lmp10/dane.lel", "r");
    read_pts_failed(in, &pts);
    make_spl(&pts, &spl);

    return 0;
}
