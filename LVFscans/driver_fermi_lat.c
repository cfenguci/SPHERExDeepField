#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include "src/spitzer_mosaic.h"

extern struct wcsprm *wcs;

int main(int argc, char *argv[])
{
  int i,j,k;
  int ii,jj;
  char (*list)[FILE_SIZE];
  int numfiles;
  int nx,ny,nz;


  nx=atoi(argv[1]);
  ny=atoi(argv[2]);
  nz=atoi(argv[3]);

  char full_mosaic_name[1024];
  strcpy(full_mosaic_name,"/home/cmb/computing/spherex/LVFscans/output/cexamples/lat_source_zmax90_gt1gev_ccube.fits");
  printf("%d\t%d\t%s\t%s\n",nx,ny,full_mosaic_name);


  numfiles=1;

  struct MOSAIC *ms=malloc(sizeof(struct MOSAIC)*numfiles);
  for(i=0;i<numfiles;i++)
  {
    ms[i].ra=create_grid(nx,ny);
    ms[i].dec=create_grid(nx,ny);
    ms[i].image_3d=create_3d_grid(nx,ny,nz);
    ms[i].exposure=create_grid(nx,ny);
    ms[i].pixx=create_grid(nx,ny);
    ms[i].pixy=create_grid(nx,ny);

    ms[i].nx=nx;
    ms[i].ny=ny;
  }


  double **image=create_grid(nx,ny);
  double value_header;
  for(i=0;i<numfiles;i++)
  {
    printf("%d\n",i);
    extract_wcs_coords(0,1,full_mosaic_name,&ms[i]);
    get_image_3d(full_mosaic_name,ms[i].image_3d);
    printf("Loaded full mosaic\n");
  }

  FILE *fp;
  char filename[]="test_fermi.fits";

  fp=fopen(filename,"w");
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) 
  {
    double val=ms[0].image_3d[i][j][6];

    if(!isnan(val)&&!isinf(val)) fprintf(fp,"%d\t%d\t%f\t%f\t%f\n",i+1,j+1,ms[0].ra[i][j],ms[0].dec[i][j],val);
    else fprintf(fp,"%d\t%d\t%f\t%f\t%f\n",i+1,j+1,ms[0].ra[i][j],ms[0].dec[i][j],0);
  }
  fclose(fp);
  



  return 0;
}

