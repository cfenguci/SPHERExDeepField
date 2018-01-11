#ifndef _UTIL_GRIDXY_H_
#define _UTIL_GRIDXY_H_


struct GRIDXY
{
  int nx,ny;
  double xmin,xmax;
  double ymin,ymax;
  double dx,dy;
  double **grid;
  double **grid_raw;
  double **hits;


};

double get_grid_spacing(const int num, const double xmin, const double xmax);
int insert_grid(const double x, const double y, const double v, GRIDXY *gxy);
void flush_grid(GRIDXY *gxy);

GRIDXY * NEW_GRIDXY(const int nx, const int ny, const double xmin, const double xmax, const double ymin, const double ymax);

#endif
