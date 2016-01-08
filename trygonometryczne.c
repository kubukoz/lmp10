#include "makespl.h"
#include "gaus/matrix.h"

void make_spl(points_t *pts, spline_t *spl) {
    int parts = 6;
    int x_size = pts->n;
    double *X = pts->x;
    double *Y = pts->y;


    //jednak nie będzie gradient descent
    //trzeba będzie wyzerowac pochodne i podać to do rozwiązania solverowi

}


//debugging
int main(int argc, char *argv[]) {
    points_t pts;
    spline_t spl;
    pts.n = 0;
    spl.n = 0;
    FILE *in = fopen("test/dane.1", "r");
    read_pts_failed(in, &pts);

    make_spl(&pts, &spl);

    return 0;
}
