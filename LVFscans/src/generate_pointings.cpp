#include "lvf_scans.h"

#include "global.h"


void generate_pointings_t0()
{
  double xmin=-3.5/2,xmax=3.5/2;
  double ymin=-3.5/2, ymax=3.5/2;
  int num=24;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j;
  double x,y;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T0");
  for(i=0;i<num;i++)
  {
    x=xmin;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  for(i=0;i<num;i++)
  {
    x=xmax;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  out.close();

}

void generate_pointings_t1()
{
  double xmin=-3.5/2,xmax=3.5/2;
  double ymin=-3.5, ymax=3.5;
  int num=48;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j;
  double x,y;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T1");
  for(i=0;i<num;i++)
  {
    x=xmin;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  for(i=0;i<num;i++)
  {
    x=xmax;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  out.close();

}

void generate_pointings_t2()
{
  double xmin=-3.5/2,xmax=3.5/2;
  double ymin=-7.0, ymax=7.0;
  int num=24;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j;
  double x,y;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T2");
  for(i=0;i<num;i++)
  {
    x=xmin+dx*i;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  out.close();

}

void generate_pointings_t3()
{
  double xmin=-5.25,xmax=5.25;
  double ymin=-7.0, ymax=7.0;
  int num=72;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j;
  double x,y;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T3");
  for(i=0;i<num;i++)
  {
    x=xmin+dx*i;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  for(i=0;i<num;i++)
  {
    x=xmax-dx*i;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }

  out.close();

}



void generate_pointings_jamie()
{
  double xmin=-3.5/2,xmax=3.5/2;
  double ymin=-3.5/2, ymax=3.5/2;

  int num=24;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j;
  double x,y;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_jamie");
  for(i=0;i<num;i++)
  {
    x=xmin;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  for(i=0;i<num;i++)
  {
    x=xmax;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  //out.close();

  int num_ext=24;
  double spacing=-dy*6;

  for(i=0;i<num_ext;i++)
  {
    x=0;
    y=ymin-dy*i-spacing;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  //out.close();

  for(i=0;i<num_ext;i++)
  {
    x=0;
    y=ymax+dy*i+spacing;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  out.close();
}



void generate_pointings_jamie_ext()
{
  double xmin=-3.5/2,xmax=3.5/2;
  double ymin=-3.5/2, ymax=3.5/2;

  int num=24;
  double dx=(xmax-xmin)/num;
  double dy=(ymax-ymin)/num;
  int i,j,ii;
  double x,y;

  int num_ext=24;
  double spacing=-dy*6;
  int n_repeat;
  n_repeat=2;

  double r,phi,theta;

  ofstream out;

  out.open("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T5");
  for(i=0;i<num;i++)
  {
    x=xmin;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;

  }
  for(i=0;i<num;i++)
  {
    x=xmax;
    y=ymin+dy*i;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  //out.close();


for(ii=0;ii<n_repeat;ii++)
{

  for(i=0;i<num_ext;i++)
  {
    x=0;
    y=ymin-dy*i-spacing;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
  //out.close();

  for(i=0;i<num_ext;i++)
  {
    x=0;
    y=ymax+dy*i+spacing;
    x*=PI/180.0;
    y*=PI/180.0;

    r=sqrt(x*x+y*y);
    theta=asin(r);
    phi=atan2(y,x);

    out<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<phi<<endl;
  }
}

  out.close();
}
