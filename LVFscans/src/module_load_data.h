
#ifndef _load_data_h_
#define _load_data_h_

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "string.h"

struct struct_file_dims
{
  int nrow, ncol;
  string filename;
};


struct struct_data_table
{
  int nrow,ncol;
  string filename;
  double **table;
};


struct_file_dims * get_file_dims(string fileName);
struct_data_table * load_data(string file);

#endif

