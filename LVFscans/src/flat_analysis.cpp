#include "lvf_scans.h"
#include "global.h"

void load_hits_file(string path, string filename, const int nx, const int ny, double **grid)
{
  int i,j;
  
  string fullname=path+filename;
  ifstream in;
//cout<<fullname<<endl;
  int idx,idy;
  int vhits;

  in.open(fullname.c_str());
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) 
  {
    double sum;
    in>>idx>>idy>>vhits>>sum;
    grid[idx][idy]=vhits;
//	cout<<idx<<"\t"<<idy<<"\t"<<vhits<<endl;
  }
  in.close();
  cout<<"Loaded hits file "<<fullname<<endl;
}





void get_angular_counts(string filename)
{
  int i,j;
  
  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  
  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  double r=(xmax-xmin)/2.0;
  double dr=sqrt(dx*dx+dy*dy);
  int numr=r/dr+1;

  double **pole_count=create_2d_grid<double>(numx,numy);
  double *rings=new double [numr+1];
  double *rings_count=new double [numr+1];
  
  //cout<<"mem"<<endl;


  for(i=0;i<numx;i++) for(j=0;j<numy;j++) pole_count[i][j]=0.0;
  for(i=0;i<numr;i++) 
  {	  
    rings[i]=0.0;
	rings_count[i]=0.0;
  }
  
  //cout<<1<<endl;

  load_hits_file(path,filename,numx,numy,pole_count);

  
  //get the rings
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
	y=ymin+j*dy;
	double radius=sqrt(x*x+y*y);
	int idr=radius/dr; 
        if(idr<numr)
        {
	  rings[idr]+=pole_count[i][j];
	  rings_count[idr]++;
        }
  }
  
  for(i=0;i<numr;i++) 
  {
    double v=rings_count[i];
	if(v!=0) rings[i]/=v;
	else rings[i]=0;
  }
  

  ofstream out;
  out.open((path+filename+"-rings").c_str());
  for(i=0;i<numr;i++) 
  {
	double r=i*dr;
    out<<i<<"\t"<<r<<"\t"<<rings[i]<<endl;
  }
  out.close();

  //delete_2d_grid<double>(pole_count,numx,numy);
  //delete rings;
  //delete rings_count;
}

void get_angular_counts(string path, string filename)
{
  int i,j;
  
  //string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  
  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  double r=(xmax-xmin)/2.0;
  double dr=sqrt(dx*dx+dy*dy);
  int numr=r/dr+1;

  double **pole_count=create_2d_grid<double>(numx,numy);
  double *rings=new double [numr+1];
  double *rings_count=new double [numr+1];
  
  //cout<<"mem"<<endl;


  for(i=0;i<numx;i++) for(j=0;j<numy;j++) pole_count[i][j]=0.0;
  for(i=0;i<numr;i++) 
  {	  
    rings[i]=0.0;
	rings_count[i]=0.0;
  }
  
  //cout<<1<<endl;

  load_hits_file(path,filename,numx,numy,pole_count);

  
  //get the rings
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
	y=ymin+j*dy;
	double radius=sqrt(x*x+y*y);
	int idr=radius/dr; 
        if(idr<numr)
        {
	  rings[idr]+=pole_count[i][j];
	  rings_count[idr]++;
        }
  }
  
  for(i=0;i<numr;i++) 
  {
    double v=rings_count[i];
	if(v!=0) rings[i]/=v;
	else rings[i]=0;
  }
  

  ofstream out;
  out.open((path+filename+"-rings").c_str());
  for(i=0;i<numr;i++) 
  {
	double r=i*dr;
    out<<i<<"\t"<<r<<"\t"<<rings[i]<<endl;
  }
  out.close();

  //delete_2d_grid<double>(pole_count,numx,numy);
  //delete rings;
  //delete rings_count;
}








//for systematic tests



map_pairs* get_angular_counts_sys(const double theta_a,const double theta_b,string filename, double **image, double **hits)
{
  int i,j;
  
  //string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  
  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  double r=(xmax-xmin)/2.0;
  double dr=sqrt(dx*dx+dy*dy);
  int numr=r/dr+1;


  double *rings=new double [numr+1];
  double *rings_count=new double [numr+1];
  

  for(i=0;i<numr;i++) 
  {	  
    rings[i]=0.0;
	rings_count[i]=0.0;
  }


  
  //get the rings
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
	y=ymin+j*dy;
	double radius=sqrt(x*x+y*y);
	int idr=radius/dr; 
    if(idr<numr)
    {
	  rings[idr]+=image[i][j];
	  rings_count[idr]++;
    }
  }
  
  for(i=0;i<numr;i++) 
  {
    double v=rings_count[i];
	if(v!=0) rings[i]/=v;
	else rings[i]=0;
  }

  ofstream out;

  map_pairs *mp=new map_pairs;
  mp->dark_current_template=create_2d_grid<double>(numx,numy);
  mp->taper=create_2d_grid<double>(numx,numy);
  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    mp->dark_current_template[i][j]=0.0;
    mp->taper[i][j]=0.0;
  }

  out.open((filename+"-radial-template").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
    y=ymin+j*dy;
    double radius=sqrt(x*x+y*y);
    int idr=radius/dr; 
    if(idr<numr) mp->dark_current_template[i][j]=rings[idr];   
    out<<i<<"\t"<<j<<"\t"<<mp->dark_current_template[i][j]<<"\t"<<0<<endl;
  }
  out.close();

  out.open((filename+"-radial-taper").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++)
  {
    x=xmin+i*dx;
    y=ymin+j*dy;
    double radius=sqrt(x*x+y*y);
    //int idr=radius/dr;
    double cutoff1=theta_a*PI/180.0;
    double cutoff2=theta_b*PI/180.0;
    if((radius>=cutoff1)&&(radius<=cutoff2)) mp->taper[i][j]=cos(PI/2.0*(radius-cutoff1)/((theta_b-theta_a)*PI/180.0));
    if(radius<cutoff1) mp->taper[i][j]=1.0;
    if(radius>cutoff2) mp->taper[i][j]=0.0;
    out<<i<<"\t"<<j<<"\t"<<mp->taper[i][j]<<"\t"<<0<<endl;
  }
  out.close();


  out.open((filename+"-rings").c_str());
  for(i=0;i<numr;i++) 
  {
	double r=i*dr;
    out<<i<<"\t"<<r<<"\t"<<rings[i]<<endl;
  }
  out.close();


  delete rings;
  delete rings_count;

  return mp;
}










void get_grid(const int id, Healpix_Map<double> &map, string name)
{
  int i,j;
  string file;
  //file="/home/cmb/computing/spherex/VLFscans/output/VLF_scan.fits";
  hpint64 nside=map.Nside(),idpix;
  //Healpix_Map<double> map(nside,RING,SET_NSIDE);
  //read_Healpix_map_from_fits(file,map,1,2);
  vec3 vec;

  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  double **pole_count=create_2d_grid<double>(numx,numy);
  double **pole_count0=create_2d_grid<double>(numx,numy);

  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    pole_count[i][j]=0.0;
    pole_count0[i][j]=0.0;
  }

  double theta,phi;

  ofstream out;
  string null="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%d",id);
  cout<<id<<endl;




  out.open((null+name+"_"+seg+"-scatter.txt").c_str());
  for(i=0;i<map.Npix();i++)
  {
    pix2ang_ring64(nside, i, &theta, &phi);
    theta/=deg2rad;
    phi/=deg2rad;

    //if(theta>=(180.0-POLE_BOUND))
    if(theta<=POLE_BOUND)
    {
      vec=map.pix2vec(i);
      int idx=(vec.x-xmin)/dx;
      int idy=(vec.y-ymin)/dy;
      pole_count[idx][idy]+=map[i];
      pole_count0[idx][idy]++;

      double x,y,v;
      x=vec.x;//sin(theta)*cos(phi);
      y=vec.y;//sin(theta)*sin(phi);
      v=map[i];
      out<<x<<"\t"<<y<<"\t"<<v<<endl;
    }

  }
  out.close();

  out.open((null+name+"_"+seg+".txt").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    if(pole_count0[i][j]<=0) out<<i<<"\t"<<j<<"\t"<<0<<endl;
    else out<<i<<"\t"<<j<<"\t"<<pole_count[i][j]<<endl;
  }
  out.close();

  delete_2d_grid<double>(pole_count,numx,numy);
  delete_2d_grid<double>(pole_count0,numx,numy);
}

void get_original_scan()
{
    get_grid(0,_map_GAL_, "grid_original");
}

