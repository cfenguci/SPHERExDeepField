#include <stdio.h>

#include <fitsio.h>

#include <wcslib.h>
#include <getwcstab.h>

#include "spitzer_mosaic.h"


double ** create_grid(const int nx, const int ny)
{
  int i,j;
  double **arr;
  arr=malloc(nx*sizeof(double *));
  for(i=0;i<nx;i++) arr[i]=malloc(ny*sizeof(double));
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) arr[i][j]=0.0;
  return arr;
}

