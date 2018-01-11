#include <stdio.h>

#include <fitsio.h>

#include <wcslib.h>
#include <getwcstab.h>

#include "spitzer_mosaic.h"


struct wcsprm *wcs;


struct pair get_bounds(const int numfiles, struct MOSAIC *ms, TYPE_MOSAIC tm1, TYPE_MOSAIC tm2)
{
  int i,j;
  int ii,jj;
  double max_v1=-1e100;
  double min_v1=1e100;
  double max_v2=-1e100;
  double min_v2=1e100;
  int nx,ny;
  double v1,v2;

  nx=ms[0].nx;
  ny=ms[0].ny;
  //printf("%d\t%d\n",nx,ny);
  for(i=0;i<numfiles;i++)
  {
    for(ii=0;ii<nx;ii++) for(jj=0;jj<ny;jj++)
    {
      switch(tm1)
      {
        case ra:
          v1=ms[i].ra[ii][jj];
        break;
        case dec:
          v1=ms[i].dec[ii][jj];
        break;
        case image:
          v1=ms[i].image[ii][jj];
        break;
        case exposure:
          v1=ms[i].exposure[ii][jj];
        break;
        case pixx:
          v1=ms[i].pixx[ii][jj];
        break;
        case pixy:
          v1=ms[i].pixy[ii][jj];
        break;
      }

      switch(tm2)
      {
        case ra:
          v2=ms[i].ra[ii][jj];
        break;
        case dec:
          v2=ms[i].dec[ii][jj];
        break;
        case image:
          v2=ms[i].image[ii][jj];
        break;
        case exposure:
          v2=ms[i].exposure[ii][jj];
        break;
        case pixx:
          v2=ms[i].pixx[ii][jj];
        break;
        case pixy:
          v2=ms[i].pixy[ii][jj];
        break;
      }

//printf("%f\t%f\n",v1,v2);
//cout<<ii<<"\t"<<jj<<"\t"<<vra<<"\t"<<vdec<<"\t"<<MOSAIC[i].image[ii][jj]<<"\t"<<MOSAIC[i].expo[ii][jj]<<endl;
      if(v1>=max_v1) max_v1=v1;
      if(v1<=min_v1) min_v1=v1;
      if(v2>=max_v2) max_v2=v2;
      if(v2<=min_v2) min_v2=v2;
    }
  }
  struct pair rect;
  rect.x_min=min_v1;
  rect.x_max=max_v1;
  rect.y_min=min_v2;
  rect.y_max=max_v2;

  printf("%f\t%f\n",rect.x_min,rect.x_max);
  printf("%f\t%f\n",rect.y_min,rect.y_max);
  return rect;
}

double ** create_grid(const int nx, const int ny)
{
  int i,j;
  double **arr;
  arr=malloc(nx*sizeof(double *));
  for(i=0;i<nx;i++) arr[i]=malloc(ny*sizeof(double));
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) arr[i][j]=0.0;
  return arr;
}


double ***create_3d_grid(const int nx, const int ny, const int nz)
{
  int i,j,k;
  double *** mat=(double ***)calloc(nx,sizeof(double **));
  for(i=0;i<nx;i++)
  {
    mat[i]=(double **)calloc(ny,sizeof(double *));
    for(j=0;j<ny;j++) mat[i][j]=(double *)calloc(nz,sizeof(double));
  }
  for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++) mat[i][j][k]=0.0;
  return mat;
}


int read_cbcd_list(char *filename, char (*list)[FILE_SIZE])
{
  FILE *fp;
//  const char filename[] = "/data-3/cmb/SpizterData/code_spitzer_ciber/cbcd_list_epoch1";
  int cnt=0;
  fp= fopen (filename, "r");
  if ( fp != NULL )
  {
    char line[FILE_SIZE]; /* or other suitable maximum line size */
    while (fgets(line,sizeof(line),fp) != NULL ) /* read a line */
    {
      //fputs(line,stdout); /* write the line */
      memcpy(list[cnt++],line, FILE_SIZE);
    }
    fclose(fp);
   }
   else
   {
     perror( filename ); /* why didn't the file open? */
   }

  return cnt;
}

int extract_wcs_coords(const int dohdr, const int dopixel, char *file,struct MOSAIC *ms)
{
//printf("%s\n",file);
  char *header, *hptr;
//  int  dohdr = 0, dopixel = 0;
  int doworld = 0;
  int  i, j, nkeyrec, nreject, nwcs, stat[NWCSFIX], status = 0;
  double imgcrd[2], phi, pixcrd[2], theta, world[2];
  fitsfile *fptr;
//  struct wcsprm *wcs;
  int nx,ny;

  nx=ms[0].nx;
  ny=ms[0].ny;
//  nx=atoi(argv[2]);
//  ny=atoi(argv[3]);

  /* Parse options. */
/*
  for (i = 1; i < argc && argv[i][0] == '-'; i++) {
    if (!argv[i][1]) break;

    switch (argv[i][1]) {
    case 'h':
      dohdr = 1;
      break;
    case 'p':
      dopixel = 1;
      break;
    case 'w':
      doworld = 1;
      break;
    default:
      fprintf(stderr, "Usage: twcshdr [-h | -p | -w] <file>\n");
      return 1;
    }
  }
*/
/*
  if (i != (argc-1)) {
    fprintf(stderr, "Usage: twcshdr [-h | -p | -w] <file>\n");
    return 1;
  }
*/

//printf("%s\n",argv[0]);
//printf("%s\n",argv[1]);
//printf("%s\n",argv[2]);
//printf("%s\n",argv[3]);

//dopixel=1;


char filein[96];
memcpy(filein,file,95);
//printf("%s\n",filein); //for(i=0;i<95;i++) printf("%c\t",filein[i]);
  /* Open the FITS test file and read the primary header. */
//char *filexx="/data-3/cmb/SpizterData/nep/epoch1/r49723392/ch1/bcd/SPITZER_I1_49723392_0080_0000_1_cbcd.fits";
//printf("%d\n",sizeof(*filexx));
//printf("%d\n",sizeof(*file));
filein[94]='\0';
//for(i=0;i<95;i++) printf("%c\t%c\n",filein[i],file[i]);

  fits_open_file(&fptr,filein, READONLY, &status);
  if ((status = fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status))) {
    fits_report_error(stderr, status);
    return 1;
  }
//printf("File openned\n");

  /*-----------------------------------------------------------------------*/
  /* Basic steps required to interpret a FITS WCS header, including -TAB.  */
  /*-----------------------------------------------------------------------*/

  /* Parse the primary header of the FITS file. */
  if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs,
                       &wcs))) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
  }


  /* Read coordinate arrays from the binary table extension. */
  if ((status = fits_read_wcstab(fptr, wcs->nwtb, (wtbarr *)wcs->wtb,
                                 &status))) {
    fits_report_error(stderr, status);
    return 1;
  }

  /* Translate non-standard WCS keyvalues. */
  if ((status = wcsfix(7, 0, wcs, stat))) {
    for (i = 0; i < NWCSFIX; i++) {
      if (stat[i] > 0) {
        fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
                wcsfix_errmsg[stat[i]]);
      }
    }

    return 1;
  }

  /*-----------------------------------------------------------------------*/
  /* The wcsprm struct is now ready for use.                               */
  /*-----------------------------------------------------------------------*/

  /* Finished with the FITS file. */
  fits_close_file(fptr, &status);
  free(header);

  /* Initialize the wcsprm struct, also taking control of memory allocated by
   * fits_read_wcstab(). */
  if ((status = wcsset(wcs))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
    return 1;
  }

  if (dohdr) 
  {
    if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) 
    {
      return 1;
    }

    hptr = header;
    printf("\n\n");
    for (i = 0; i < nkeyrec; i++, hptr += 80) 
    {
      printf("%.80s\n", hptr);
    }
    free(header);
  } 
  else if (dopixel) 
  {
//printf("Do pixel %d\t%d\n",nx,ny);
    for(i=0;i<nx;i++) for(j=0;j<ny;j++)
    {
      //printf("Enter pixel coordinates: ");
      //if (scanf("%lf%*[ ,]%lf", pixcrd, pixcrd+1) != wcs->naxis) break;
      pixcrd[0]=i;
      pixcrd[1]=j;
      status = wcsp2s(wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
      (*ms).ra[i][j]=world[0];
      (*ms).dec[i][j]=world[1];
      //printf("  (%20.15f, %20.15f) ->\n  (%20.15f, %20.15f)\n\n", pixcrd[0], pixcrd[1], world[0], world[1]);
//printf("%20.15f\t%20.15f\t%20.15f\t%20.15f\n",pixcrd[0], pixcrd[1], world[0], world[1]);
    }

  } 
  else if (doworld) 
  {
    while (1) 
    {
      printf("Enter world coordinates: ");
      if (scanf("%lf%*[ ,]%lf", world, world+1) != wcs->naxis) break;
      status = wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
      //printf("  (%20.15f, %20.15f) ->\n  (%20.15f, %20.15f)\n\n", world[0], world[1], pixcrd[0], pixcrd[1]);
    }

  } 
  else 
  {
    /* Print the struct. */
    if ((status = wcsprt(wcs))) 
    {
      return 1;
    }
  }

  /* Clean up. */
  //status = wcsvfree(&nwcs, &wcs);
  return 0;
}



int generate_pixel_coords(struct MOSAIC ms)
{
//printf("%s\n",file);
  char *header, *hptr;
  int  dohdr = 0, dopixel = 0, doworld = 0;
  int  i, j, nkeyrec, nreject, nwcs, stat[NWCSFIX], status = 0;
  double imgcrd[2], phi, pixcrd[2], theta, world[2];
  fitsfile *fptr;
//  struct wcsprm *wcs;
  int nx,ny;
  nx=ms.nx;
  ny=ms.ny;
//  nx=atoi(argv[2]);
//  ny=atoi(argv[3]);

  for(i=0;i<nx;i++) for(j=0;j<ny;j++)
  {
    //printf("Enter world coordinates: ");
    //if (scanf("%lf%*[ ,]%lf", world, world+1) != wcs->naxis) break;
    world[0]=ms.ra[i][j];
    world[1]=ms.dec[i][j];
    status = wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
//    printf("  (%20.15f, %20.15f) ->\n  (%20.15f, %20.15f)\n\n", world[0], world[1], pixcrd[0], pixcrd[1]);
    ms.pixx[i][j]=pixcrd[0];
    ms.pixy[i][j]=pixcrd[1];
  }


  return 0;
}



