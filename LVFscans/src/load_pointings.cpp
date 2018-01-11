
#include "global.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "string.h"


double** load_pointings(string file)
{
  int len=144;
  int i,j;
  double **arr_pointings=create_2d_grid<double>(len,2);

  double x,y,theta,phi;

  ifstream in;

  in.open(file.c_str());
  for(i=0;i<len;i++)
  {
    in>>x>>y>>theta>>phi;
    arr_pointings[i][0]=theta;
    arr_pointings[i][1]=phi;
  }
  in.close();
  cout<<"Loaded "<<file<<endl;

  return arr_pointings;
}
