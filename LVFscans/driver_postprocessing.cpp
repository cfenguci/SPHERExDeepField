#include "src/lvf_scans.h"
#include "src/global.h"




int main(int argc, char *argv[])
{
  hpint64 nside=512;
  Healpix_Map<double> _map_GAL_(nside,RING,SET_NSIDE);

  //load_fits_map_cxx("/home/cfeng/computing/spherex/LVFscans/output/wmap_band_imap_r9_9yr_W_v5.fits",scanned_map);
  postprocessing("NULL");
  return 1;
}
