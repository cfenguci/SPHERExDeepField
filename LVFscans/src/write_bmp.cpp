#include "lvf_scans.h"

void WriteBMP(const char *fname, int w,int h,unsigned char *img)
{
    FILE *f = fopen(fname,"wb");

    unsigned char bfh[54] = {0x42, 0x4d,
    /* bfSize [2]*/ 54, 0, 0, 0, /**/
    /* reserved [6]*/ 0, 0, 0, 0, /**/
    /* biOffBits [10]*/ 54, 0, 0, 0, /**/
    /* biSize [14]*/ 40, 0, 0, 0, /**/
    /* width [18]*/ 0, 0, 0, 0, /**/
    /* height [22]*/ 0, 0, 0, 0, /**/
    /* planes [26]*/ 1, 0, /**/
    /* bitcount [28]*/ 24, 0,/**/
    /* compression [30]*/ 0, 0, 0, 0, /**/
    /* size image [34]*/ 0, 0, 0, 0, /**/
    /* xpermeter [38]*/ 0, 0, 0, 0, /**/
    /* ypermeter [42]*/ 0, 0, 0, 0, /**/
    /* clrused [46]*/ 0, 0, 0, 0, /**/
    /* clrimportant [50]*/ 0, 0, 0, 0 /**/};
    int realw = w * 3, rem = w % 4, isz = (realw + rem) * h, fsz = isz + 54;
    //bfh.bfSize = fsz;
    bfh[2] = (fsz & 0xFF); bfh[3] = (fsz >> 8) & 0xFF; bfh[4] = (fsz >> 16) & 0xFF; bfh[5] = (fsz >> 24) & 0xFF;
    //bfh.biSize = isz
    bfh[34] = (isz & 0xFF); bfh[35] = (isz >> 8) & 0xFF; bfh[36] = (isz >> 16) & 0xFF; bfh[37] = (isz >> 24) & 0xFF;
    //bfh.biWidth = w;
    bfh[18] = (w & 0xFF); bfh[19] = (w >> 8) & 0xFF; bfh[20] = (w >> 16) & 0xFF; bfh[21] = (w >> 24) & 0xFF;
    //bfh.biHeight = h;
    bfh[22] = (h & 0xFF); bfh[23] = (h >> 8) & 0xFF; bfh[24] = (h >> 16) & 0xFF; bfh[25] = (h >> 24) & 0xFF;

    // xpels/ypels
    // bfh[38] = 19; bfh[39] = 11;
    // bfh[42] = 19; bfh[43] = 11;

    fwrite((void*)bfh, 54, 1, f);

    unsigned char* bstr = new unsigned char[realw], *remstr = 0; 
    if(rem != 0) { remstr = new unsigned char[rem]; memset(remstr,0,rem); }
/*
    for(int j = h-1 ; j > -1 ; j--){
            for(int i = 0 ; i < w ; i++)
                    for(int k = 0 ; k < 3 ; k++) { bstr[i*3+k] = img[(j*realw+i*3)+(2-k)]; }
            fwrite(bstr,realw,1,f); if (rem != 0) { fwrite(remstr,rem,1,f); }
    }
*/
    //for(int j = h-1 ; j > -1 ; j--)
    for(int j =0 ; j<=h-1; j++)
    {
      for(int i = 0 ; i < w ; i++)
      for(int k = 0 ; k < 3 ; k++) 
      { 
        bstr[i*3+k] = img[(j*realw+i*3)+(2-k)]; 
      }
      fwrite(bstr,realw,1,f); if (rem != 0) { fwrite(remstr,rem,1,f); }
    }


    delete [] bstr; 
    if(remstr) delete [] remstr;

    fclose(f);
}


