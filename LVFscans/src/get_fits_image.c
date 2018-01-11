#include <string.h>
#include <stdio.h>
#include "fitsio.h"


#include "spitzer_mosaic.h"

int get_image(char *file, double **image)
{
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis, ii, anynul;
    long naxes[2] = {1,1}, fpixel[2] = {1,1};
    double *pixels;
    char format[20], hdformat[20];

char filein[96];
memcpy(filein,file,95);
//printf("%s\n",filein); //for(i=0;i<95;i++) printf("%c\t",filein[i]);
  /* Open the FITS test file and read the primary header. */
//char *filexx="/data-3/cmb/SpizterData/nep/epoch1/r49723392/ch1/bcd/SPITZER_I1_49723392_0080_0000_1_cbcd.fits";
//printf("%d\n",sizeof(*filexx));
//printf("%d\n",sizeof(*file));
filein[94]='\0';


    if (!fits_open_file(&fptr, filein, READONLY, &status))
    {
        //int hdunum,hdutype,status;
        //hdunum=idim;
        //if(fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) print("Can not move to HDU %d\n",idim);

        if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) )
        {
          if (naxis > 2 || naxis == 0)
             printf("Error: only 1D or 2D images are supported\n");
          else
          {
            /* get memory for 1 row */
            pixels = (double *) malloc(naxes[0] * sizeof(double));

            if (pixels == NULL) {
                printf("Memory allocation error\n");
                return(1);
            }

            if (bitpix > 0) {  /* set the default output format string */
               strcpy(hdformat, " %7d");
               strcpy(format,   " %7.0f");
            } else {
               strcpy(hdformat, " %15d");
               strcpy(format,   " %15.5f");
            }

              /* terminate header line */

            /* loop over all the rows in the image, top to bottom */
            for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
            {
               if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL,
                    pixels, NULL, &status) )  /* read row of pixels */
                  break;  /* jump out of loop on error */

               //printf(" %4d ",fpixel[1]);  /* print row number */
               for (ii = 0; ii < naxes[0]; ii++)
			   //{
                  //printf(format, pixels[ii]);
                  image[ii][fpixel[1]-1]=pixels[ii];
			   //}
			  /* print each value  */
               //printf("\n");                    /* terminate line */
            }
            free(pixels);
          }
        }
        fits_close_file(fptr, &status);
    } 

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
