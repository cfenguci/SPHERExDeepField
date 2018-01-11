
#include "lvf_scans.h"


void equ2gal(const double ra, const double dec, double &l, double &b)
{
        double a_NGP, delta_NGP, l_ascend, x, y;

        a_NGP=192.8594813;
        delta_NGP=27.1282511;
        l_ascend=33;

        a_NGP*=(PI/180.0);
        delta_NGP*=(PI/180.0);
        l_ascend*=(PI/180.0);




        b=sin(dec)*sin(delta_NGP)+cos(dec)*cos(delta_NGP)*cos(ra-a_NGP);
        y=sin(dec)*cos(delta_NGP)-cos(dec)*sin(delta_NGP)*cos(ra-a_NGP);
        x=cos(dec)*sin(ra-a_NGP);

        b=asin(b);
        l=atan2(y, x);

        l+=l_ascend;

}

void gal2equ(const double l, const double b, double &ra, double &dec)
{
        double a_NGP, delta_NGP, l_ascend, x, y, z;

        a_NGP=192.8594813;
        delta_NGP=27.1282511;
        l_ascend=33;

        a_NGP*=(PI/180.0);
        delta_NGP*=(PI/180.0);
        l_ascend*=(PI/180.0);

        z=sin(b)*sin(delta_NGP)+cos(b)*cos(delta_NGP)*sin(l - l_ascend);
        dec=asin(z);

        y=cos(b)*cos(l - l_ascend);
        x=sin(b)*cos(delta_NGP) - cos(b)*sin(delta_NGP)*sin(l - l_ascend);
        ra=atan2(y,x)+a_NGP;

}













