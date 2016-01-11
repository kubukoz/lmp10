#include "makespl.h"
#include "gaus/matrix.h"
#include <math.h>
#include <stdlib.h>
#include "gaus/piv_ge_solver.h"

#define DEFAULT_PARTS 10;

typedef double (*FUNC_DOUBLE)(double);

int get_n(int i) {
    return (i - 1) / 2 + 1;
}

double a_func(double x) {
    return 1;
}

FUNC_DOUBLE get_func(int i) {
    if (i == 0) return a_func;
    else if (i % 2 == 1) return sin;
    else return cos;
}

double f(double *a, int n, double x) {
    FUNC_DOUBLE func;
    double value = 0;
    int i;

    for (i = 0; i < n; i++) {
        int curr_n = get_n(i);
        func = get_func(i);
        value += a[i] * func(x * curr_n);
    }
    return value;
}

double f1(double *a, int n, double x) {
    FUNC_DOUBLE func;
    double value = 0;
    int i;

    /*od i=1, bo wszystkie pochodne ze stałej są zerowe*/
    for (i = 1; i < n; i++) {
        int curr_n = get_n(i);
        /*przejscie funkcji*/
        func = get_func(i + 1);
        /*zmiana znaku dla cos*/
        int sign = i % 2 == 0 ? -1 : 1;
        value += a[i] * sign * curr_n * func(curr_n * x);
    }
    return value;
}

double f2(double *a, int n, double x) {
    FUNC_DOUBLE func;
    double value = 0;
    int i;
    for (i = 1; i < n; i++) {
        int curr_n = get_n(i);
        func = get_func(i);
        /*w drugiej pochodnej zmieniamy znaki obu funkcji, stąd -= */
        value -= a[i] * pow(curr_n, 2) * func(curr_n * x);
    }
    return value;
}

double f3(double *a, int n, double x) {
    FUNC_DOUBLE func;
    double value = 0;
    int i;
    for (i = 1; i < n; i++) {
        int curr_n = get_n(i);
        /*przejscie funkcji*/
        func = get_func(i + 1);
        /*zmiana znaku dla sin*/
        int sign = i % 2 != 0 ? -1 : 1;
        value += a[i] * sign * pow(curr_n, 3) * func(curr_n * x);
    }
    return value;
}

double *get_factors(matrix_t *mx) {
    int i;
    double *results = malloc(mx->rn * sizeof(double));
    for (i = 0; i < mx->rn; i++) {
        results[i] = get_entry_matrix(mx, i, mx->cn - 1);
    }
    return results;
}

void make_spl(points_t *pts, spline_t *spl) {
    int default_parts = DEFAULT_PARTS;

    int i, j, k;

    char *env = getenv("APPROX_BASE_SIZE");
    int parts = env == NULL ? default_parts : atoi(env);

    int x_size = pts->n;

    matrix_t *equations = make_matrix(parts, parts + 1);

    /*zmienna k odpowiada za rząd - czyli dzięki k zmienia się: POCHODNA.*/
    for (k = 0; k < parts; k++) {
        FUNC_DOUBLE derivative = get_func(k);

        /*n takie, że derivative=sin/cos(n*x)*/
        int derN = get_n(k);

        /*zmienna j odpowiada za kolumnę - czyli zmienia: FUNKCJĘ.*/
        for (j = 0; j < parts + 1; j++) {
            /*suma w kolumnie*/
            double sum = 0;
            FUNC_DOUBLE func = get_func(j);

            /*n takie, że func=sin/cos(n*x)*/
            double funcN = get_n(j);

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

    if (piv_ge_solver(equations)) {
        spl->n = 0;
        return;
    }

    double *factors = get_factors(equations);

    free(equations->e);
    free(equations);

    /*tworzymy spline'y*/
    if(alloc_spl(spl, pts->n + 1)){
        spl->n = 0;
        return;
    }

    double start_x = pts->x[0];
    double range = pts->x[x_size - 1] - start_x;
    double delta = range * 1.0 / (pts->n - 1);
    for (i = 0; i < pts->n + 1; i++) {
        double current_x = start_x + delta * i;
        spl->x[i] = current_x;
        spl->f[i] = f(factors, parts, current_x);
        spl->f1[i] = f1(factors, parts, current_x);
        spl->f2[i] = f2(factors, parts, current_x);
        spl->f3[i] = f3(factors, parts, current_x);
    }

    free(factors);
}
