#ifndef _bin_data_h_
#define _bin_data_h_

struct UNBINNED_DATA
{
int len;
double *x,*y;
};

struct BINNED_DATA
{
int len;
double bin_start;
double bin_step;
double *bin;
double *bin_mid;
double *bin_value;
};

UNBINNED_DATA * NEW_UNBINNED_DATA(const int len);
BINNED_DATA *NEW_BINNED_DATA(const int len, const double binstart, const double binstep);
void bin_data(UNBINNED_DATA *ud, BINNED_DATA *bd);




struct band_info
{
  int num_bin;
  int num_band;
  double *bin;
};



struct band_averaged
{
  int num_bin;
  int num_band;
  double *bin;
  double *band_mid;
  double *band_value;

};

band_info * NEW_band_info(const int num_bin);
band_averaged * NEW_band_averaged(const int num_bin);
band_averaged * bin_data(UNBINNED_DATA *ud, band_info *bi);


#endif
