#include "healpix_include.h"

double wave_2_freq(const double lambda);
double freq_2_wave(const double nu);

void INIT_NCP_flat();
double get_list_minimum(const int len, double *arr);
void Analyze_NCP_flat(const int);
void Analyze_NCP_image(const int day);

void Generate_NCP_noise(const int, Healpix_Map<double> &,Healpix_Map<double> &,Healpix_Map<double> &);
void RECORD_NCP_flat(const int whichband, double *FOV);
void capture_NCP_image(const int,const double vsky, double *FOV);

void OUTPUT_NCP_flat(const int);

double ** load_sphere_band_noise(string filename);

void compute_cl_from_flat(const int,const double,const double,const int,const int whichone);
