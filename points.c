#include "points.h"
#include <stdlib.h>

int
realloc_pts_failed (points_t * pts, int size)
{
  return (pts->x = realloc (pts->x, size * sizeof *pts->x)) == NULL
      || (pts->y = realloc(pts->y, size * sizeof *pts->y)) == NULL;
}

int
read_pts_failed (FILE * inf, points_t * pts)
{
  int size;
  double x, y;

  if (pts->n == 0) {
    pts->x = malloc (sizeof(double));
    if (pts->x == NULL)
      return 1;
    pts->y = malloc (sizeof(double));
    if (pts->y == NULL) {
      free (pts->x);
      return 1;
    }
    size = 1;
  }
  else
    size = pts->n + 1;

  while (fscanf (inf, "%lf %lf", &x, &y) == 2) {
    pts->x[pts->n] = x;
    pts->y[pts->n] = y;
    pts->n++;
    if (!feof (inf)) {
      if (realloc_pts_failed (pts, size + 1))
        return 1;
      else
        size++;
    }
  }

  if (pts->n != size)
    if (realloc_pts_failed (pts, pts->n))
      return 1;

  return 0;
}
