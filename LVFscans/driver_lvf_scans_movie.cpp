#include "src/lvf_scans.h"


extern GLint windW,windH;
extern int _num_orbit_;

int main(int argc, char *argv[])
{
  int i;  

  int total_num=atoi(argv[1]);
  for(i=0;i<total_num;i++)
  {
    _num_orbit_=i;

    get_position(0);
  }


  return 1;
}
