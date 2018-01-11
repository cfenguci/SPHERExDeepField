#include "lvf_scans.h"

extern hpint64 _nside_;
extern Healpix_Map<double> _map_NAIVE_;
extern Healpix_Map<double> _map_LVF_scan_;
extern Healpix_Map<double> _map_GAL_;

void make_mask(const double theta_0, Healpix_Map<double> &mask)
{
  double theta,phi;
  hpint64 nside=_nside_;
  int i;
  for(i=0;i<mask.Npix();i++)
  {
    pix2ang_ring64(nside, i, &theta, &phi);
    theta/=deg2rad;
    phi/=deg2rad;

    if(theta<theta_0) mask[i]=1;
    else mask[i]=0;
  }
}

void get_naive_map_powspec()
{
  int i,j;
  ofstream out;
  double l,fac;
  double fsky;
  hpint64 nside=_nside_;
  int nlmax,nmmax;
  nlmax=2*nside;
  nmmax=nlmax;
  arr<double> weight;
  paramfile params;
  int num_iter=0;


  string file;


  Healpix_Map<double> map_original(nside,RING,SET_NSIDE);
  Healpix_Map<double> map(nside,RING,SET_NSIDE);
  Healpix_Map<double> mask(nside,RING,SET_NSIDE);

  Alm<xcomplex<double> > alm_original(nlmax,nmmax);
  Alm<xcomplex<double> > alm(nlmax,nmmax);
  PowSpec powspec_original,powspec_scan;

  make_mask(45.0,mask);

  get_ring_weights(params,map.Nside(),weight);

  int npix=mask.Npix();

  for(i=0;i<npix;i++) 
  {
    double v0,v1;
    v0=_map_GAL_[i]*mask[i];
    v1=_map_NAIVE_[i]*mask[i];

    map[i]=v1;
    map_original[i]=v0;

    if(isnan(v0)||isinf(v0)||(fabs(v0)>100)) 
    {
      map_original[i]=0.0;
      mask[i]=0.0;
    }
    if(isnan(v1)||isinf(v1)||(fabs(v1)>100)) 
    {
//cout<<i<<"\t"<<map[i]<<"\t";
      map[i]=0.0;
      mask[i]=0.0;
//cout<<map[i]<<endl;
    }

  }

//for(i=0;i<npix;i++) cout<<i<<"\t"<<map[i]<<endl;



  double sum=0.0;
  for(i=0;i<npix;i++) sum+=mask[i];
  fsky=sum/npix;
  cout<<"Fsky="<<fsky<<endl;

  map2alm_iter(map_original,alm_original,num_iter,weight);
  map2alm_iter(map,alm,num_iter,weight);

  extract_crosspowspec(alm_original, alm_original, powspec_original);
  extract_crosspowspec(alm, alm, powspec_scan);
  
  file="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/spherex_naive_map_cl";


  out.open(file.c_str());
  for(i=0;i<nlmax;i++)
  {
    l=i;
    fac=l*(l+1.0)/2.0/PI;

    powspec_original.tt(i)/=fsky;
    powspec_scan.tt(i)/=fsky;
    double cl1,cl2;
    cl1=powspec_original.tt(i);
    cl2=powspec_scan.tt(i);

    out<<l<<"\t"<<fac*cl1<<"\t"<<cl1<<"\t"<<fac*cl2<<"\t"<<cl2<<endl;
//    out<<l<<"\t"<<powspec_original.tt(i)<<"\t"<<powspec_scan.tt(i)<<endl;
  }
  out.close();


}

