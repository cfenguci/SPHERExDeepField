

#include "src/lvf_scans.h"


#include "src/global.h"
#include "src/flat_ncp.h"

/*
double wave_2_freq(const double wave)
{
//  double LIGHT_SPEED=299792458
  double freq=LIGHT_SPEED/wave*1e6;
  return freq;
}
*/


int main(int argc, char *argv[])
{
  int whichgrid=atoi(argv[1]);
  int id_dprj=atoi(argv[2]);
  double theta_a=atof(argv[3]);
  double theta_b=atof(argv[4]);
  compute_cl_from_flat(id_dprj,theta_a,theta_b,whichgrid,360);
  return 1;
}
