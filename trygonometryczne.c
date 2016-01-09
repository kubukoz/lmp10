#include "makespl.h"
#include "gaus/matrix.h"
#include <math.h>
#include <stdlib.h>
#include "gaus/piv_ge_solver.h"

#define DEFAULT_PARTS 5;

double a_func(double x) {
    return 1;
}

double *get_factors(matrix_t *);

void make_spl(points_t *pts, spline_t *spl) {
    int default_parts = DEFAULT_PARTS;

    int i, j, k;

    char *env = getenv("APPROX_BASE_SIZE");
    int parts = env == NULL ? default_parts : atoi(env);

    int x_size = pts->n;

    matrix_t *equations = make_matrix(parts, parts + 1);

    /*zmienna k odpowiada za rząd - czyli dzięki k zmienia się: POCHODNA.*/
    for (k = 0; k < parts; k++) {
        double (*derivative)(double);

        /*n takie, że derivative=sin/cos(n*x)*/
        int derN = (k-1) / 2 + 1;

        if(k == 0) derivative = &a_func;
        else if (k % 2 == 1) derivative = &sin;
//        if(k % 2 == 0 ) derivative = &sin;
        else derivative = &cos;

        /*zmienna j odpowiada za kolumnę - czyli zmienia: FUNKCJĘ.*/
        for (j = 0; j < parts + 1; j++) {
            /*suma w kolumnie*/
            double sum = 0;
            double (*func)(double);

            /*n takie, że func=sin/cos(n*x)*/
            double funcN = (j-1) / 2 + 1;

            if(j == 0) func = &a_func;
            else if (j % 2 == 1) func = &sin;
//            if(j % 2 == 0) func = &sin;
            else func = &cos;

            if (j < parts) {
                /*dodajemy pochodną razy wartość funkcji bazowej dla danego x*/
                for (i = 0; i < x_size; i++) {
                    double x = pts->x[i];

                    sum += func(x * funcN)
                           * derivative(x * derN);
                }
            }
            else {
                /*dodajemy pochodną dla x razy wartość y które znamy*/
                for (i = 0; i < x_size; i++) {
                    double x = pts->x[i];
                    double y = pts->y[i];

                    sum += y * derivative(x * derN);
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
