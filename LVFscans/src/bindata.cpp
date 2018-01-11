#include "bindata.h"

UNBINNED_DATA * NEW_UNBINNED_DATA(const int len)
{
  UNBINNED_DATA *p=new UNBINNED_DATA;
  p->len=len;
  p->x=new double [len];
  p->y=new double [len];

  return p;
}

/*
int len;
double bin_start;
double bin_step;
double *bin;
double *bin_mid;
double *bin_value;
*/
BINNED_DATA *NEW_BINNED_DATA(const int len, const double binstart, const double binstep)
{
  BINNED_DATA *p=new BINNED_DATA;
  p->len=len;
  p->bin_start=binstart;
  p->bin_step=binstep;
  p->bin=new double [len];
  p->bin_mid=new double [len];
  p->bin_value=new double [len];

  return p;
}


void bin_data(UNBINNED_DATA *ud, BINNED_DATA *bd)
{
  int i,j;

  for(i=0;i<bd->len;i++) bd->bin[i]=bd->bin_start+bd->bin_step*i;
  for(i=0;i<bd->len-1;i++) bd->bin_mid[i]=(bd->bin[i]+bd->bin[i+1])/2.0;
  double inte;
  for(i=0;i<bd->len-1;i++)
  {
    inte=0.0;
    for(j=0;j<ud->len;j++) if((ud->x[j]-bd->bin[i])*(ud->x[j]-bd->bin[i+1])<=0) inte+=ud->y[j];
    bd->bin_value[i]=inte/(bd->bin[i+1]-bd->bin[i]);

  }
}



band_info * NEW_band_info(const int num_bin)
{
  band_info *bi=new band_info;
  bi->num_bin=num_bin;
  bi->num_band=num_bin-1;
  bi->bin=new double [num_bin];

  return bi;
}


band_averaged * NEW_band_averaged(const int num_bin)
{
  band_averaged *ba=new band_averaged;
  ba->num_bin=num_bin;
  ba->num_band=num_bin-1;
  ba->bin=new double [num_bin];
  ba->band_value=new double [num_bin-1];
  ba->band_mid=new double [num_bin-1];

  return ba;
}



band_averaged * bin_data(UNBINNED_DATA *ud, band_info *bi)
{
  int i,j;
  int nbin=bi->num_bin;
  int nband=bi->num_band;

  band_averaged *p=NEW_band_averaged(nbin);

  double inte;
  double x,x1,current_x;
  for(i=0;i<nband;i++)
  {
    inte=0.0;
    x=bi->bin[i];
    x1=bi->bin[i+1];
    for(j=0;j<ud->len;j++)
    {
      current_x=ud->x[j];
      if((current_x>=x)&&(current_x<x1)) inte+=ud->y[j];
    }
    p->band_value[i]=inte/(x1-x);
    p->band_mid[i]=(x+x1)/2.0;
  }

  return p;
}


