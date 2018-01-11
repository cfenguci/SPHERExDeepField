#include <stdio.h>

#include <fitsio.h>

#include <wcslib.h>
#include <getwcstab.h>


#define NX 256//8503
#define NY 256//8521

#define FILE_SIZE 1024

//#define FULL_MOSAIC_MAP

struct MOSAIC
{
  int nx,ny;
  double **ra;//[NX][NY];
  double **dec;//[NX][NY];
  double **image;//[NX][NY];
  double ***image_3d;

  double **exposure;//[NX][NY];

  double **pixx;//[NX][NY];
  double **pixy;//[NX][NY];

};

typedef enum{ra,dec,image,exposure,pixx,pixy} TYPE_MOSAIC;

struct pair
{
  double x_min,x_max;
  double y_min,y_max;
};

//struct wcsprm *wcs;

int get_image(char *file, double **image);
int get_image_3d(char *file, double ***image);

int get_header(char *file,char *name_header_in, double *);

void pixel2radec();

struct pair get_bounds(const int numfiles, struct MOSAIC *ms, TYPE_MOSAIC tm1, TYPE_MOSAIC tm2);
double ** create_grid(const int nx, const int ny);
double ***create_3d_grid(const int nx, const int ny, const int nz);
int read_cbcd_list(char *,char (*list)[FILE_SIZE]);
int extract_wcs_coords(const int dohdr, const int dopixel, char *file,struct MOSAIC *ms);
int generate_pixel_coords(struct MOSAIC ms);
