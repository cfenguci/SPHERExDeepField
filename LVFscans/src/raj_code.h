#ifndef _raj_code_h_
#define _raj_code_h_

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


struct_file_dims *get_file_dims(string file);
struct_data_table *load_data_table(string file);

void get_grid(string, string);


#endif
