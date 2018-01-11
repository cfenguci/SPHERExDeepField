
#include "lvf_scans.h"
#include "global.h"
#include "raj_code.h"
#include "bindata.h"

double wave_2_freq(const double lambda)
{
  return LIGHT_SPEED/lambda*1e6;//[Hz]
}
double freq_2_wave(const double nu)
{
return LIGHT_SPEED/nu*1e6;//um
}



void INIT_NCP_flat()
{
  int i;
  int nx,ny;
  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;

  nx=801;//801;
  ny=801;//801;

  GDXY=new GRIDXY *[48];
  for(i=0;i<48;i++) GDXY[i]=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);

  image_ncp_s=new GRIDXY;
  image_ncp_s=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);

  image_ncp_n=new GRIDXY;
  image_ncp_n=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);

  image_ncp_zodi=new GRIDXY;
  image_ncp_zodi=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);


  image_ncp_dark=new GRIDXY;
  image_ncp_dark=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);


  image_ncp_thermo=new GRIDXY;
  image_ncp_thermo=NEW_GRIDXY(nx,ny,xmin,xmax,ymin,ymax);


}

void RECORD_NCP_flat(const int whichband, double *FOV)
{
  int i;
  double x=FOV[0];
  double y=FOV[1];
  double z=FOV[2];
  double v=1.0;

  insert_grid(x,y,v,GDXY[whichband]);
}

void capture_NCP_image(const int whichgrid,const double vsky, double *FOV)
{
  int i;
  double x=FOV[0];
  double y=FOV[1];
  double z=FOV[2];

  switch(whichgrid)
  {
    case 0:
      insert_grid(x,y,vsky,image_ncp_s);
    break;

    case 1:
      insert_grid(x,y,vsky,image_ncp_n);
    break;


    case 2:
      insert_grid(x,y,vsky,image_ncp_zodi);
    break;

    case 3:
      insert_grid(x,y,vsky,image_ncp_dark);
    break;

    case 4:
      insert_grid(x,y,vsky,image_ncp_thermo);
    break;


  }
}

double get_list_minimum(const int len, double *arr)
{
  int i;
  double min=1e100;
  double v;
  for(i=0;i<len;i++) 
  {
    v=arr[i]; 
    if(v<=min) min=v;
  }
  return min;
}

double get_list_sum(const int len, double *arr)
{
  double sum=0.0;
  int i;
  for(i=0;i<len;i++) sum+=arr[i];
  return sum;
}

void Analyze_NCP_flat(const int day)
{
  GRIDXY *gd=GDXY[0];
  int nx=gd->nx;
  int ny=gd->ny;
  int i,j,idband;

  int len=48;
  double *arr=new double [len];
  double least,sum;

  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%03d",day);

  string file;
  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/LeastHits_";
  file=path+seg;
  ofstream out;

  //double **hits=create_2d_grid<double>(nx,ny);

  out.open(file.c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++)
  {
    for(idband=0;idband<len;idband++) arr[idband]=GDXY[idband]->grid_raw[i][j];
    least=get_list_minimum(len,arr);
    sum=get_list_sum(len,arr);
    //hits[i][j]=least;
    out<<i<<"\t"<<j<<"\t"<<least<<"\t"<<sum<<endl;
  }
  out.close();
  delete arr;
/*
  string frame="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/fov_frame_";
  for(idband=0;idband<len;idband++)
  {
    char seg_fov[1024];
    memset(seg_fov,0,1024);
    sprintf(seg_fov,"%02d",idband);
    string frame_name=frame+seg+"-id_"+seg_fov;
    out.open(frame_name.c_str());
    for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<GDXY[idband]->grid_raw[i][j]<<endl;
    out.close();
  }
*/
  //return hits;

}


void Analyze_NCP_image(const int day)
{
  
  int nx=image_ncp_s->nx;
  int ny=image_ncp_s->ny;
  int i,j,idband;



  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%03d",day);

  string file;
  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/Flat_image_";
  file=path;
  ofstream out;

  flush_grid(image_ncp_s);
  flush_grid(image_ncp_n);
  flush_grid(image_ncp_zodi);
  flush_grid(image_ncp_dark);
  flush_grid(image_ncp_thermo);

  out.open((file+"s-"+seg).c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<image_ncp_s->grid[i][j]<<"\t"<<image_ncp_s->hits[i][j]<<endl;
  out.close();
  out.open((file+"n-"+seg).c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<image_ncp_n->grid[i][j]<<"\t"<<image_ncp_n->hits[i][j]<<endl;
  out.close();
  out.open((file+"zodi-"+seg).c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<image_ncp_zodi->grid[i][j]<<"\t"<<image_ncp_zodi->hits[i][j]<<endl;
  out.close();
  out.open((file+"dark-"+seg).c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<image_ncp_dark->grid[i][j]<<"\t"<<image_ncp_dark->hits[i][j]<<endl;
  out.close();
  out.open((file+"thermo-"+seg).c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) out<<i<<"\t"<<j<<"\t"<<image_ncp_thermo->grid[i][j]<<"\t"<<image_ncp_thermo->hits[i][j]<<endl;
  out.close();

}



void Generate_NCP_noise(const int idband, Healpix_Map<double> &map, Healpix_Map<double> &mask, Healpix_Map<double> &Nobs)
{

  string filename,path;

  int i,j;

  int len=48;

  filename="/home/cfeng/computing/spherex/LVFscans/output/Spherex_band_noise";
  double **band_noise=load_sphere_band_noise(filename);

  //for(i=0;i<48;i++) cout<<band_noise[i][1]<<"\t"<<band_noise[i][2]<<endl;

  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;


  double **pole_count=create_2d_grid<double>(numx,numy);
  path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  filename="LeastHits_360";
  load_hits_file(path,filename,numx,numy,pole_count);



  double *arr=new double [len];
  double least;
  char seg[1024];
  double eta,noise;
  double *noise_level=new double [len];
  for(i=0;i<len;i++) 
  {
    double r=i+1.0;
    noise_level[i]=band_noise[i][2];
  }

  string file;
  path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/Map_band_noise_BANDID-";

  hpint64 nside=map.Nside(),idpix;
  for(i=0;i<mask.Npix();i++) 
  {
    mask[i]=0.0;
    map[i]=0.0;
  }
  cout<<"Clean mask and map"<<endl;

  double theta,phi,r;


  ofstream out;


  double t1=4,t2=5.5;
  double **mask_edge=create_2d_grid<double>(numx,numy);
  for(i=0;i<numx;i++) for(j=0;j<numy;j++) mask_edge[i][j]=0.0;
  //out.open((path+"-mask").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
    y=ymin+j*dy;
    double radius=sqrt(x*x+y*y)/deg2rad;
    //int idr=radius/dr;
    if(radius<t1) mask_edge[i][j]=1.0;
    if((radius>=t1)&&(radius<=t2)) mask_edge[i][j]=cos(PI/2.0*(radius-t1)/(t2-t1));
    if(radius>t2) mask_edge[i][j]=0.0;
    //if(idr<numr) mp->dark_current_template[i][j]=rings[idr];
    //out<<i<<"\t"<<j<<"\t"<<mask_edge[i][j]<<"\t"<<0<<endl;

    //mask_edge[i][j]=1.0;
  }
  //out.close();




  cout<<"Creating noise Band="<<idband<<"\t Sigma="<<noise_level[idband]<<" [kJy/sr]"<<endl;




  //for(idband=0;idband<1;idband++)
  {
    memset(seg,0,1024);
    sprintf(seg,"%02d",idband);
    file=path+seg;
    //out.open(file.c_str());
    for(i=0;i<numx;i++) for(j=0;j<numy;j++)
    {
      x=xmin+i*dx;
      y=ymin+j*dy;

      r=sqrt(x*x+y*y);
      theta=asin(r);
      phi=atan2(y,x);

      ang2pix_ring64(nside,theta,phi,&idpix);

      eta=_rng_.rand_gauss();
      double nhits=pole_count[i][j];

      if(nhits>0) 
      {
        //nhits=330;
        noise=eta*noise_level[idband]/sqrt(nhits);
        //mask[idpix]=mask_edge[i][j];
        Nobs[idpix]=noise_level[idband]*noise_level[idband]/nhits;
      }
      else 
      {
        noise=0;//eta*noise_level[idband]; 
        //mask[idpix]=0;
        Nobs[idpix]=0;
      }

      map[idpix]=noise;
      mask[idpix]=mask_edge[i][j];

      //out<<i<<"\t"<<j<<"\t"<<x<<"\t"<<y<<"\t"<<noise<<endl;

    }
    //out.close();
  

  }


  delete arr;


}

void OUTPUT_NCP_flat(const int day)
{
  int i,j,idband;
  char seg[1024];

  char sday[1024];
  memset(sday,0,1024);
  sprintf(sday,"%03d",day);

  string file;
  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/BandHits_Day-";
//  file=path+sday+"_Band-";
  ofstream out;

  for(idband=0;idband<48;idband++)
  {
    memset(seg,0,1024);
    sprintf(seg,"%d",idband);
    file=path+sday+"_Band-"+seg;

    GRIDXY *gd=GDXY[idband];
    
    out.open(file.c_str());
    for(i=0;i<gd->nx;i++) for(j=0;j<gd->ny;j++) out<<i<<"\t"<<j<<"\t"<<gd->grid_raw[i][j]<<endl;
    out.close();
  }
}

double ** load_sphere_band_noise(string filename)
{
  ifstream in;
  double **table_sphere_band_noise=create_2d_grid<double>(48,4);
  int i,j;
  double wave,sigma,sigma2;
  int id;
  in.open(filename.c_str());
  for(i=0;i<48;i++) 
  {
    in>>id>>wave>>sigma>>sigma2;

    table_sphere_band_noise[id][0]=id;
    table_sphere_band_noise[id][1]=wave;
    table_sphere_band_noise[id][2]=sigma;
    table_sphere_band_noise[id][3]=sigma2;
  }
  in.close();

  return table_sphere_band_noise;
}




void compute_cl_from_flat(const int id_dprj,const double theta_a,const double theta_b,const int whichgrid,const int whichone)
{

  string filename,path;
  //string *name={"s","n","zodi","dark","thermo"};
  vector<string> name(5);
  name[0]="s";
  name[1]="n";
  name[2]="zodi";
  name[3]="dark";
  name[4]="thermo";

  int i,j;

  char seg[1024];
  sprintf(seg,"%03d",whichone);

  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  hpint64 nside=1024;
  Healpix_Map<double> map(nside,RING,SET_NSIDE);
  Healpix_Map<double> mask(nside,RING,SET_NSIDE);
  hpint64 idpix;


  double **pole_count=create_2d_grid<double>(numx,numy);
  path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  filename=path+"Flat_image_"+name[whichgrid]+"-"+seg;
  struct_data_table *dt=load_data_table(filename);
  cout<<"Loaded "<<filename<<endl;



  int nx=801;
  double **mat=create_2d_grid<double>(nx,nx);
  double **mat_hits=create_2d_grid<double>(nx,nx);

  for(i=0;i<dt->nrow;i++)
  {
    int idx=dt->table[i][0];
    int idy=dt->table[i][1];
    mat[idx][idy]=dt->table[i][2];
    mat_hits[idx][idy]=dt->table[i][3];
  }


  map_pairs *mp;
  if(id_dprj) 
  {
    cout<<"Get dark current radial distrubution"<<endl;
    mp=get_angular_counts_sys(theta_a,theta_b,filename,mat,mat_hits);
  
    for(i=0;i<nx;i++) for(j=0;j<nx;j++) 
    {
      mat[i][j]-=mp->dark_current_template[i][j];
      mat_hits[i][j]=mp->taper[i][j];
    }

  }

  for(i=0;i<mask.Npix();i++) 
  {
    mask[i]=0.0;
    map[i]=0.0;
  }
  cout<<"Clean mask and map"<<endl;

  double theta,phi,r;

  ofstream out;

  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
    y=ymin+j*dy;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    ang2pix_ring64(nside,theta,phi,&idpix);
      
    double vsky=mat[i][j];
    double vhits;
    vhits=mat_hits[i][j];

    if(vhits>0)    
    {
      map[idpix]=vsky;
      if(id_dprj) mask[idpix]=vhits;
      else mask[idpix]=1.0;
     
    }
  }



  int npix=map.Npix();
  double pixel=4.0*PI/npix;
//get the DC
  double vDC=0.0;
  double npix_eff=0.0;
  for(i=0;i<npix;i++)  
  {
    vDC+=map[i];
    if(mask[i]>0) npix_eff+=1.0;
  }
  vDC/=npix_eff;
  cout<<"DC level="<<vDC<<endl;

  for(i=0;i<npix;i++) map[i]-=vDC;

  



  double fsky=0.0;
  double sum=0.0;
  for(i=0;i<npix;i++) sum+=mask[i];
  fsky=sum/npix;
  cout<<"Sky fraction is "<<fsky<<"\t or \t"<<4*PI*(180.0/PI)*(180.0/PI)*fsky<<" deg^2"<<endl;

  arr<double> weights;

  paramfile params;
  int num_iter=0;
  int nlmax=2048;
  int nmmax=nlmax;
  Alm<xcomplex<double> > alm(nlmax,nmmax);
  PowSpec powspec(1,nlmax);
  get_ring_weights(params,map.Nside(),weights);


  for(i=0;i<npix;i++) map[i]*=mask[i];



  map2alm_iter(map,alm,num_iter,weights);
  extract_crosspowspec(alm, alm, powspec);




  filename=path+"cl_Flat_image_"+name[whichgrid]+"-"+seg;
  out.open(filename.c_str());

  for(i=0;i<nlmax;i++)
  {
    double cl=powspec.tt(i);
    double l,fac;
    l=i;
    fac=l*(l+1.0)/2.0/PI;

    cl/=fsky;
    out<<l<<"\t"<<fac*cl<<"\t"<<cl<<endl;
  }
  out.close();


/*
  int len_bin=20;//32;
  double bin_start=100;
  double bin_step=100;

  UNBINNED_DATA *ud=NEW_UNBINNED_DATA(nlmax);
  BINNED_DATA *bd=NEW_BINNED_DATA(len_bin, bin_start, bin_step);

  for(i=0;i<ud->len;i++)
  {
    double fac,l;
    l=i;
    ud->x[i]=l;
    ud->y[i]=powspec.tt(i)/fsky; 
  }
  bin_data(ud,bd);
  filename=path+"bin_Flat_image_"+name[whichgrid]+"-"+seg;
  out.open(filename.c_str());
  for(i=0;i<bd->len;i++)
  {
    double l=bd->bin_mid[i];
    double fac=l*(l+1.0)/2.0/PI;
    out<<l<<"\t"<<fac*bd->bin_value[i]<<"\t"<<bd->bin_value[i]<<endl;
  }
  out.close();
*/

  int len_bin=10;//32;
  double bin_start=100;
  double bin_step=100;
  UNBINNED_DATA *ud=NEW_UNBINNED_DATA(nlmax);

//init a log binning
  band_info *bi=NEW_band_info(len_bin+1);
  double l0=100.0;
  double delta=(log10(2200)-log10(l0))/len_bin;
  for(i=0;i<=len_bin;i++) bi->bin[i]=l0*pow(10.0,i*delta);
  //bi->bin[len_bin]=2000.0;

  for(i=0;i<ud->len;i++)
  {
    double fac,l;
    l=i;
    ud->x[i]=l;
    ud->y[i]=powspec.tt(i)/fsky; 
  }

  band_averaged *p=bin_data(ud,bi);

  filename=path+"bin_Flat_image_"+name[whichgrid]+"-"+seg;
  out.open(filename.c_str());
  for(i=0;i<p->num_band;i++)
  {
    double l=p->band_mid[i];
    double fac=l*(l+1.0)/2.0/PI;
    double bp=p->band_value[i];

    out<<l<<"\t"<<fac*bp<<"\t"<<bp<<endl;
  }


}



