
#include "lvf_scans.h"
#include "util_gridxy.h"

double get_grid_spacing(const int num, const double xmin, const double xmax)
{
  return (xmax-xmin)/num;
}

int insert_grid(const double x, const double y, const double v, GRIDXY *gxy)
{
  double xmin=gxy->xmin;
  double xmax=gxy->xmax;
  double ymin=gxy->ymin;
  double ymax=gxy->ymax;
  double dx=gxy->dx;
  double dy=gxy->dy;
  int flag=0;

  int idx,idy;
  idx=(x-xmin)/dx;
  idy=(y-ymin)/dy;

  if( (idx>=0)&&(idx<gxy->nx)&&(idy>=0)&&(idy<gxy->ny) ) 
  {
    flag=1;
    gxy->grid_raw[idx][idy]+=v;
    gxy->hits[idx][idy]++;

  }
  else flag=0;


}

void flush_grid(GRIDXY *gxy)
{
  int nx=gxy->nx;
  int ny=gxy->ny;
  int i,j;

  for(i=0;i<nx;i++) for(j=0;j<ny;j++)
  {
    double v=gxy->hits[i][j];
    if(v!=0) gxy->grid[i][j]=gxy->grid_raw[i][j]/v;
    else gxy->grid[i][j]=0.0;
  }

}

GRIDXY * NEW_GRIDXY(const int nx, const int ny, const double xmin, const double xmax, const double ymin, const double ymax)
{
  GRIDXY *p=new GRIDXY;

  p->nx=nx;
  p->ny=ny;
  p->xmin=xmin;
  p->xmax=xmax;
  p->ymin=ymin;
  p->ymax=ymax;

  p->dx=get_grid_spacing(nx,xmin,xmax);
  p->dy=get_grid_spacing(ny,ymin,ymax);

  p->grid=create_2d_grid<double>(nx,ny);
  p->grid_raw=create_2d_grid<double>(nx,ny);
  p->hits=create_2d_grid<double>(nx,ny);

  int i,j;
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) 
  {
    p->grid[i][j]=0.0;
    p->grid_raw[i][j]=0.0;
    p->hits[i][j]=0.0;
  }




  return p;

}
