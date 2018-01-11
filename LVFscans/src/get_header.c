#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include "spitzer_mosaic.h"


int get_header(char *file,char *name_header_in, double *value_header)
{
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD], newcard[FLEN_CARD];
    char oldvalue[FLEN_VALUE], comment[FLEN_COMMENT];
    int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
    int iomode, keytype;
int i;
char filein[96];
char NAME_header[8];
memcpy(filein,file,95);
//printf("%s\n",filein); //for(i=0;i<95;i++) printf("%c\t",filein[i]);

  /* Open the FITS test file and read the primary header. */
//char *filexx="/data-3/cmb/SpizterData/nep/epoch1/r49723392/ch1/bcd/SPITZER_I1_49723392_0080_0000_1_cbcd.fits";
//printf("%d\n",sizeof(*filexx));
//printf("%d\n",sizeof(*file));
filein[94]='\0';

memcpy(NAME_header,name_header_in,8);
NAME_header[7]='\0';
//for(i=0;i<8;i++) printf("%c\n",NAME_header[i]);
/*
    if (argc == 3) 
      iomode = READONLY;
    else if (argc == 4)
      iomode = READWRITE;
    else {
      printf("Usage:  modhead filename[ext] keyword newvalue\n");
      printf("\n");
      printf("Write or modify the value of a header keyword.\n");
      printf("If 'newvalue' is not specified then just print \n");
      printf("the current value. \n");
      printf("\n");
      printf("Examples: \n");
      printf("  modhead file.fits dec      - list the DEC keyword \n");
      printf("  modhead file.fits dec 30.0 - set DEC = 30.0 \n");
      return(0);
    }
*/
    if (!fits_open_file(&fptr, filein, iomode, &status))
    {
      if (fits_read_card(fptr,NAME_header, card, &status))
      {
        printf("Keyword does not exist\n");
        card[0] = '\0';
        comment[0] = '\0';
        status = 0;  /* reset status after error */
      }
      //else
       // printf("%s\n",card);
//for(int k;k<FLEN_CARD;k++)
//printf("%d\t%c\n",k,card[k]);
      char * pch;
        //printf ("Splitting string \"%s\" into tokens:\n",str);
//double value_header;
int cnt=0;
//value_header=atof(pch);
//printf("%.6f\n",value_header);
        pch = strtok (card,"=/");
cnt++;
        while (pch != NULL)
        {
            //printf ("%d\t%s\n",cnt,pch);
if(cnt==2)
{
*value_header=atof(pch);
//printf("%.6f\n",value_header);
}
pch = strtok (NULL, "=/");
cnt++;
        }
      //if (argc == 4)  /* write or overwrite the keyword */
      {
          /* check if this is a protected keyword that must not be changed */
          if (*card && fits_get_keyclass(card) == TYP_STRUC_KEY)
          {
            printf("Protected keyword cannot be modified.\n");
          }
          else
          {
            /* get the comment string */
            if (*card)fits_parse_value(card, oldvalue, comment, &status);

            /* construct template for new keyword */
            //strcpy(newcard, argv[2]);     /* copy keyword name */
            //strcat(newcard, " = ");       /* '=' value delimiter */
            //strcat(newcard, argv[3]);     /* new value */
            //if (*comment) {
            //  strcat(newcard, " / ");  /* comment delimiter */
            //  strcat(newcard, comment);     /* append the comment */
            //}

            /* reformat the keyword string to conform to FITS rules */
            fits_parse_template(newcard, card, &keytype, &status);

            /* overwrite the keyword with the new value */
            //fits_update_card(fptr, argv[2], card, &status);

            //printf("Keyword has been changed to:\n");
            //printf("%s\n",card);
          }
      }  /* if argc == 4 */
      fits_close_file(fptr, &status);
    }    /* open_file */

    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);
    return(status);
}

