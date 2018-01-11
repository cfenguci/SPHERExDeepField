#ifndef _systematics_h_
#define _systematics_h_

#include "raj_code.h"
void prepare_systematics();
double get_stamp_pixel(const double x, const double y, double **p);

double ** bin_field(double **p);
#endif
