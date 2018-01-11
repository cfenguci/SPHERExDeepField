#include "src/lvf_scans.h"
#include "src/global.h"
#include "src/raj_code.h"



int main(int argc, char *argv[])
{
  string file=argv[1];
//  get_grid("/home/cfeng/computing/spherex/LVFscans/output/data_points_least_deep_fig2_91.0");
  string path="/home/cfeng/computing/spherex/LVFscans/output/raj/";
//  string output="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  get_grid(path,file);
  get_angular_counts(path,file+"-grid");
  return 1;
}
