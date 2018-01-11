#include "src/lvf_scans.h"


#include "src/global.h"
#include "src/flat_ncp.h"

/*
double wave_2_freq(const double wave)
{
//  double LIGHT_SPEED=299792458
  double freq=LIGHT_SPEED/wave*1e6;
  return freq;
}
*/


int main(int argc, char *argv[])
{
  int i,idband;
  ofstream out, report, out1;
  hpint64 nside=1024;
  hpint64 idpix;
  Healpix_Map<double> map(nside,RING,SET_NSIDE),mask(nside,RING,SET_NSIDE);
  Healpix_Map<double> Nobs(nside,RING,SET_NSIDE);

  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";


  int len=48;

  string filename="/home/cfeng/computing/spherex/LVFscans/output/Spherex_band_noise_MEV";
  double **band_noise=load_sphere_band_noise(filename);
  for(i=0;i<len;i++) cout<<i<<"\t"<<band_noise[i][1]<<endl;

  report.open((path+"stat_band_noise").c_str());
  //cout<<wave_2_freq(1)<<endl;
for(idband=0;idband<48;idband++)
{
  Generate_NCP_noise(idband,map,mask,Nobs);

  int npix=map.Npix();
  double pixel=4.0*PI/npix;

  double fsky=0.0;
  double sum=0.0;
  for(i=0;i<npix;i++) sum+=mask[i];
  fsky=sum/npix;
  cout<<"Sky fraction is "<<fsky<<"\t or \t"<<4*PI*(180.0/PI)*(180.0/PI)*fsky<<" deg^2"<<endl;

  arr<double> weights;

  paramfile params;
  int num_iter=0;
  int nlmax=1024;
  int nmmax=nlmax;
  Alm<xcomplex<double> > alm(nlmax,nmmax);
  PowSpec powspec(1,nlmax);
  get_ring_weights(params,map.Nside(),weights);


  for(i=0;i<npix;i++) map[i]*=mask[i];



  double factor=0.0;
  for(i=0;i<npix;i++) if(mask[i]) factor+=mask[i]*Nobs[i];
  double theory=pixel*pixel/4.0/PI/fsky*factor;

  report<<idband<<"\t"<<theory<<"\t"<<sqrt(theory)<<"\t"<<sqrt(theory/pixel)<<endl;

  map2alm_iter(map,alm,num_iter,weights);
  extract_crosspowspec(alm, alm, powspec);

  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%02d",idband);

  //string outfile=path+"map_noise_BAND-"+seg+".fits";
  //write_Healpix_map_to_fits(outfile,map,planckType<double>());
  //cout<<"Band noise map saved "<<outfile<<endl;

  string file=path+"CL_BandNoise-"+seg;
  string file1=path+"CL_BandNoise-nW-"+seg;
  out.open(file.c_str());
  out1.open(file1.c_str());
  for(i=0;i<nlmax;i++)
  {
    double nl=powspec.tt(i);
    double wave=band_noise[idband][1];
    double conversion=1e17/wave_2_freq(wave)/1e3;
    double l,fac;
    l=i;
    fac=l*(l+1.0)/2.0/PI;

    nl/=fsky;
    out<<i<<"\t"<<nl<<"\t"<<theory<<endl;
    out1<<i<<"\t"<<nl/conversion/conversion<<"\t"<<theory/conversion/conversion<<"\t"<<wave<<"\t"<<fac*nl/conversion/conversion<<endl;
    cout<<idband<<"\t"<<i<<"\t"<<nl<<"\t"<<theory<<endl;
  }
  out.close();
  out1.close();
  cout<<"Band "<<idband<<" is done"<<endl; 

}

report.close();



  return 1;
}
