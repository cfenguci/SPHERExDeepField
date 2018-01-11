#include <stdio.h>

#include <fitsio.h>

#include <wcslib.h>
#include <getwcstab.h>

#include "spitzer_mosaic.h"

void pixel2radec()
{


const double tol = 1.0e-10;


/* In real life these would be encoded as FITS header keyrecords. */
const int NAXIS = 2;
const double CRPIX[2] =  {513.0000000,513.0000000};
const double PC[2][2] = {{0.00048747,-0.0019},
                         {0.0019,0.00048747}};
const double CRVAL[2] = {270.6741,66.3415};

const double CDELT[2] =  {0.0, 0.0};

char CTYPE[4][9] = {"WAVE-F2W", "XLAT-BON", "TIME-LOG", "XLON-BON"};

const double LONPOLE  = 180.0;
const double LATPOLE  = 0;
const double RESTFRQ  =   0.0;
const double RESTWAV  =   0.0;

int NPV = 2;
struct pvcard PV[2]; 
int i,j;
  struct wcsprm *wcss,*wcs1;
double *pcij;

char *header, *hptr;
  int  nkeyrec, nreject, nwcs, stat[NWCSFIX], status = 0;
  double imgcrd[2], phi, pixcrd[2], theta, world[2];

fitsfile *fptr;
//char *full_mosaic_name="/data-3/cmb/SpizterData/nep/epoch1/r49723392/ch1/bcd/SPITZER_I1_49723392_0080_0000_1_cbcd.fits";//"/data-3/cmb/SpizterData/nep/maps/1200mas/ch1_e1.skymap.00.fits";
char *full_mosaic_name="/data-3/cmb/SpizterData/nep/maps/ch1_half1.skymap.00.fits";
  fits_open_file(&fptr,full_mosaic_name, READONLY, &status);
  if ((status = fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status))) {
    fits_report_error(stderr, status);
    return 1;
  }


//  if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs,
//                       &wcss))) {
//    fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
//  }
//for(i=0;i<200;i++) printf("%c\n",header[i]);

int cnt0=0;
/*
char *pch;
pch = strtok (header," ");
cnt0++;
while (pch != NULL)
{

printf ("%d\t%s\n",cnt0,pch);
pch = strtok (NULL," ");
cnt0++;
}
*/

/*
char cc='0';
for(i=0;i<4960;i++)
{
printf("%d\t%c\n",i,header[i]);
cnt0++;
}
*/

header[266]='1';
header[267]='0';
header[268]='2';
header[269]='4';

header[346]='1';
header[347]='0';
header[348]='2';
header[349]='4';

header[893]='0';
header[894]='.';
header[895]='0';
header[896]='0';
header[897]='0';
header[898]='4';
header[899]='8';
header[900]='7';
header[901]='4';
header[902]='7';
header[903]='0';
header[904]='0';
header[905]='0';
header[906]='0';
header[907]='0';
header[908]='0';
header[909]='0';


header[973]='0';
header[974]='.';
header[975]='0';
header[976]='0';
header[977]='1';
header[978]='9';
header[979]='0';
header[980]='0';
header[981]='0';
header[982]='0';
header[983]='0';
header[984]='0';
header[985]='0';
header[986]='0';
header[987]='0';
header[988]='0';
header[989]='0';


header[1052]='-';
header[1053]='0';
header[1054]='.';
header[1055]='0';
header[1056]='0';
header[1057]='1';
header[1058]='9';
header[1059]='0';
header[1060]='0';
header[1061]='0';
header[1062]='0';
header[1063]='0';
header[1064]='0';
header[1065]='0';
header[1066]='0';
header[1067]='0';
header[1068]='0';
header[1069]='0';


header[1133]='0';
header[1134]='.';
header[1135]='0';
header[1136]='0';
header[1137]='0';
header[1138]='4';
header[1139]='8';
header[1140]='7';
header[1141]='4';
header[1142]='7';
header[1143]='0';
header[1144]='0';
header[1145]='0';
header[1146]='0';
header[1147]='0';
header[1148]='0';
header[1149]='0';

//crpix1
header[1217]='5';
header[1218]='1';
header[1219]='2';
header[1220]='.';
header[1221]='0';
header[1222]='0';
header[1223]='0';
header[1224]='0';
header[1225]='0';
header[1226]='0';
header[1227]='0';
header[1228]='0';
header[1229]='0';


//crpix2
header[1297]='5';
header[1298]='1';
header[1299]='2';
header[1300]='.';
header[1301]='0';
header[1302]='0';
header[1303]='0';
header[1304]='0';
header[1305]='0';
header[1306]='0';
header[1307]='0';
header[1308]='0';
header[1309]='0';

//crval1
header[1377]='2';
header[1378]='7';
header[1379]='0';
header[1380]='.';
header[1381]='6';
header[1382]='7';
header[1383]='4';
header[1384]='1';
header[1385]='0';
header[1386]='0';
header[1387]='0';
header[1388]='0';
header[1389]='0';

//crval2
header[1457]='6';
header[1458]='6';
header[1459]='.';
header[1460]='3';
header[1461]='4';
header[1462]='1';
header[1463]='5';
header[1464]='0';
header[1465]='0';
header[1466]='0';
header[1467]='0';
header[1468]='0';
header[1469]='0';


header[2109]='0';
header[2749]='0';
header[3389]='0';
header[4189]='0';

/*
char cc='0';
for(i=0;i<4960;i++)
{
printf("%d\t%c\n",i,header[i]);
cnt0++;
}
*/

//header[893]='+';
//header[973]='+';
//header[1052]='+';
/*
char *pch;
pch = strtok (header," ");
cnt0++;
while (pch != NULL)
{

printf ("%d\t%s\n",cnt0,pch);
pch = strtok (NULL," ");
cnt0++;
}
*/

  if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs,
                       &wcss))) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
  }

/*

header[893]='';
header[973]=''; 
header[1052]='';
char cc='0';
for(i=0;i<4960;i++)
{
printf("%d\t%c\n",i,header[i]);
cnt0++;
}
*/

/*
  if ((status = fits_read_wcstab(fptr, wcs->nwtb, (wtbarr *)wcs->wtb,
                                 &status))) {
    fits_report_error(stderr, status);
    return 1;
  }

  if ((status = wcsfix(7, 0, wcs, stat))) {
    for (i = 0; i < NWCSFIX; i++) {
      if (stat[i] > 0) {
        fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
                wcsfix_errmsg[stat[i]]);
      }
    }

    return 1;
  }
*/
//  fits_close_file(fptr, &status);
//  free(header);

  if ((status = wcsset(wcss))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
    return 1;
  }
//printf("%f\n",wcs->latpole);


//wcsprt(wcs);
  //int alloc=1;
  //int naxis=NAXIS;
//printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  
#ifdef _reuse_
for (j = 0; j < NAXIS; j++) {
    wcss->crpix[j] = CRPIX[j];
  }

int cnt=0;
  for (i = 0; i < NAXIS; i++) 
    for (j = 0; j < NAXIS; j++) 
{
*(wcss->pc+cnt)=0;
*(wcss->cd+cnt)=0;
cnt++;
}
//printf("%f\t%f\n",wcs->pc[0],wcs->cd[0]);
//printf("%f\t%f\n",wcs->pc[1],wcs->cd[1]);
//printf("%f\t%f\n",wcs->pc[2],wcs->cd[2]);
//printf("%f\t%f\n",wcs->pc[3],wcs->cd[3]);


  pcij = wcss->pc;
cnt=0;
  for (i = 0; i < NAXIS; i++) {
    for (j = 0; j < NAXIS; j++) {

*(wcss->pc+cnt)=PC[i][j];
*(wcss->cd+cnt)=PC[i][j];
//printf("%f\t%f\n",(*pcij),(*(wcs->pc+cnt)));
cnt++;
//pcij++;
//wcs->pc++;
    }
  }
  for (i = 0; i < NAXIS; i++) {
    wcss->cdelt[i] = CDELT[i];
  }

//  for (i = 0; i < NAXIS; i++) {
//    strcpy(wcs->ctype[i], &CTYPE[i][0]);
//  }
//printf("%s\t%s\n",wcs->ctype[0], wcs->ctype[1]);
//strcpy(wcs->ctype[0], "XLAT-BON");
//strcpy(wcs->ctype[1], "XLON-BON");


  for (i = 0; i < NAXIS; i++) {
    wcss->crval[i] = CRVAL[i];
  }

  //wcs->lonpole = wcs->crval[0];//LONPOLE;
  //wcs->latpole =  wcs->crval[1];;//LATPOLE;

  wcss->restfrq = RESTFRQ;
  wcss->restwav = RESTWAV;

wcss->altlin=2;
//wcs->lonpole =180;
//wcs->latpole=0;
(*wcss).cel.theta0=180;
(*wcss).cel.phi0=0;
wcss->cel.ref[2]=180;
wcss->cel.ref[3]=0;

/*
  wcs->npv = NPV;
  for (i = 0; i < NPV; i++) {
    wcs->pv[i] = PV[i];
  }
*/
wcss->flag=0;
  if ((status = wcsset(wcss))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
    return 1;
  }
#endif
//wcsset(wcss);
//wcs->cdelt[0]=0;
//wcs->cdelt[1]=0;
//wcs->flag=0;
//wcsset(wcs);
//wcs->types[0]=2200;
//wcs->types[1]=2201;
//wcs->crder[0]=0;
//wcs->crder[1]=0;
//wcsini(1,2,wcs);
wcsprt(wcss);

//printf("%d\n",wcs->wtb->ndim);

//r(i=0;i<1024;i++) for(j=0;j<1024;j++)
i=0;j=0;
{
pixcrd[0]=i;
pixcrd[1]=j;
status = wcsp2s(wcss, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
//if((i==512)&&(j==512))
printf("%d\t%d\t%f\t%f\t%f\t%f\n",i,j,world[0],world[1],imgcrd[0],imgcrd[1]);
}
//printf("%f\t%f\n",(wcs->lin).crpix[0],(wcs->lin).crpix[1]);
//world[0]=273.1668;
//world[1]=65.2167;
//status = wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
//printf("%f\t%f\t%f\t%f\n",pixcrd[0],pixcrd[1],world[0],world[1]);
}
