#include "makespl.h"
#include "gaus/matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gaus/piv_ge_solver.h"

#define PART_COUNT 1;

double diff_sin(double x, int n) {
    return 1 * n * cos(n * x);
}

double diff_cos(double x, int n) {
    return -1 * n * sin(n * x);
}

double *get_factors(matrix_t *);

void make_spl(points_t *pts, spline_t *spl) {
    int i, j, k;
    int parts = PART_COUNT;
    int x_size = pts->n;
    matrix_t *equations = make_matrix(parts, parts + 1);

    //rows
    for (k = 0; k < parts; k++) {
        //parts - cols
        for (j = 0; j < parts; j++) {
            double sum = 0;
            double (*func)(double, int);

            //x
            for (i = 0; i < x_size; i++) {
                int n = j / 2 + 1;
                if (j % 2 == 0) func = &diff_sin;
                else func = &diff_cos;

                sum += pts->x[i] * pts->x[i];
//                sum += func(pts->x[j], n) * pts->x[i];
            }
            put_entry_matrix(equations, k, j, sum);
        }

        //this part is ok
        double sumY = 0;
        for (i = 0; i < x_size; i++) {
            sumY += pts->y[i] * pts->x[i];
        }
        put_entry_matrix(equations, k, parts, sumY);
    }

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
    FILE *in = fopen("dane.lel", "r");
    read_pts_failed(in, &pts);
    make_spl(&pts, &spl);


    return 0;
}
