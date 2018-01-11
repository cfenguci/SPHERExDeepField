#ifndef flat_analysis_h
#define flat_analysis_h

#include "string.h"
#include "healpix_include.h"


struct map_pairs
{
double **dark_current_template;
double **taper;
};


void load_hits_file(string path, string filename, const int nx, const int ny, double **grid);
void get_angular_counts(string filename);
void get_angular_counts(string path, string filename);

//void load_sys_image_file(string path, string filename, const int nx, const int ny, double **grid);
//void get_angular_counts_sys(string path, string filename);
map_pairs* get_angular_counts_sys(const double,const double,string filename, double **image, double **hits);

void get_grid(const int id, Healpix_Map<double> &map,string name);
void get_original_scan();

#endif
