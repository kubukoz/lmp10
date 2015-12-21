#include "makespl.h"
#include "gaus/matrix.h"

void make_spl(points_t *pts, spline_t *spl) {
    int parts = 6;
    int x_size = pts->n;
    double *X = pts->x;
    double *Y = pts->y;


    //tu będzie gradient descent
    //tak z 1000 iteracji powinno na początek wystarczyć
    //wyniki funkcji błędu doda się do jakiegoś wektorka wraz z użytą theta
    //weźmie się theta dla której błąd był najmniejszy
    //jak już wyjdzie theta (parametry), to obliczamy dla 10 podpunktów
    //wartość funkcji i 3 jej pochodnych

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
