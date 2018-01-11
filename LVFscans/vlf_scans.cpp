#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "math.h"

#define PI 3.1415926535
#define deg2rad (PI/180.0)
#define MATTYPE_ZERO 0.00

template<typename MATTYPE>
MATTYPE **create_2d_grid(const int nx, const int ny)
{
  int i,j;
  MATTYPE **mat=(MATTYPE **)calloc(nx,sizeof(MATTYPE *));
  for(i=0;i<nx;i++) mat[i]=(MATTYPE *)calloc(ny,sizeof(MATTYPE));
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) mat[i][j]=(MATTYPE) MATTYPE_ZERO;
  return mat;
}

template<typename MATTYPE>
MATTYPE ***create_3d_grid(const int nx, const int ny, const int nz)
{
  int i,j,k;
  MATTYPE *** mat=(MATTYPE ***)calloc(nx,sizeof(MATTYPE **));
  for(i=0;i<nx;i++)
  {
    mat[i]=(MATTYPE **)calloc(ny,sizeof(MATTYPE *));
    for(j=0;j<ny;j++) mat[i][j]=(MATTYPE *)calloc(nz,sizeof(MATTYPE));
  }
  for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++) mat[i][j][k]=(MATTYPE) MATTYPE_ZERO;
  return mat;

}



//theta/phi deggree
int sphere_normal_vec(const double theta, const double phi, double *x)
{
  double theta_rad,phi_rad;

  theta_rad=theta*deg2rad;
  phi_rad=phi*deg2rad;

  double ct,st,cp,sp;
  ct=cos(theta_rad);
  st=sin(theta_rad);
  cp=cos(phi_rad);
  sp=sin(phi_rad);

  x[0]=-st*st*cp;
  x[1]=st*st*sp;
  x[2]=st*ct;

  double modu=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  x[0]/=modu;
  x[1]/=modu;
  x[2]/=modu;

  int flag=0;
  if(modu==0) 
  {
    flag=1;
    cout<<"Sigularity at: ("<<theta<<","<<phi<<")"<<endl;
  }
  else flag=0;
  return flag;
}

void generate_pointing(const double theta_0, const double phi_0,const double omega_theta, const double omega_phi,const double t,
double *x,double *xn)
{
  double theta,phi;
  theta=theta_0+omega_theta*t;
  phi=phi_0+omega_phi*t;
//  if(theta>=180.0) theta-=180.0;
//  if(phi>=360.0) phi-=360.0;

  //double xn[3];
  sphere_normal_vec(theta,phi,xn);
 // cout<<"normal: ("<<theta<<"\t"<<phi<<")->("<<xn[0]<<"\t"<<xn[1]<<"\t"<<xn[2]<<")"<<endl;



  x[0]=sin(deg2rad*theta)*cos(deg2rad*phi);
  x[1]=sin(deg2rad*theta)*sin(deg2rad*phi);
  x[2]=cos(deg2rad*theta);

  //return x;
} 


void vecA_cross_vecB(double *A, double *B, double *C)
{
  int i;
  C[0]=A[1]*B[2]-A[2]*B[1];
  C[1]=A[2]*B[0]-A[0]*B[2];
  C[2]=A[0]*B[1]-A[1]*B[0];
}

void AdotB(double **A,double **B, double **C, const int dimA,const int dim, const int dimB)
{
  int i,j,k;
  for(i=0;i<dimA;i++) for(k=0;k<dimB;k++) C[i][k]=0.0;

  for(i=0;i<dimA;i++) for(k=0;k<dimB;k++)
  for(j=0;j<dim;j++) C[i][k]+=A[i][j]*B[j][k];
}

int main(int argc, char *argv[])
{
  int i,j,k;
  int ii,jj;
  double t,tmin,tmax,dt;
  int num;
  double x1[3],x2[3],xdiff[3];
  double xn[3],xlong[3];

  double rot[3][3];

  tmin=0.0;
  tmax=24.0;
  dt=0.5;
  num=(tmax-tmin)/dt;


  double **Amat=create_2d_grid<double>(3,3);
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);



  for(i=100;i<101;i++)
  {
    t=tmin+dt*i;

    generate_pointing(0,0,180.0/48.0,1.0/48.0,t,x1,xn);
    generate_pointing(0,0,180.0/48.0,1.0/48.0,t+dt,x2,xn);

    for(ii=0;ii<3;ii++) xdiff[ii]=x2[ii]-x1[ii];
    double modu=sqrt(xdiff[0]*xdiff[0]+xdiff[1]*xdiff[1]+xdiff[2]*xdiff[2]);
    for(ii=0;ii<3;ii++) xdiff[ii]/=modu;

    vecA_cross_vecB(xdiff,xn,xlong);

    //cout<<t<<"\t"<<xdiff[0]<<"\t"<<xdiff[1]<<"\t"<<xdiff[2]<<endl;
    //cout<<t<<"\t"<<-xlong[0]<<"\t"<<-xlong[1]<<"\t"<<-xlong[2]<<endl;
    //cout<<t<<"\t"<<xn[0]<<"\t"<<xn[1]<<"\t"<<xn[2]<<endl;

    rot[0][0]=xdiff[0];
    rot[1][0]=xdiff[1];
    rot[2][0]=xdiff[2];

    rot[0][1]=-xlong[0];
    rot[1][1]=-xlong[1];
    rot[2][1]=-xlong[2];

    rot[0][2]=xn[0];
    rot[1][2]=xn[1];
    rot[2][2]=xn[2];
/*
    cout<<t<<endl;
    for(int ii=0;ii<3;ii++) 
    {
      for(int jj=0;jj<3;jj++) cout<<rot[ii][jj]<<"\t"; 
      cout<<endl;
    }
*/
    double w=7.0/48.0*deg2rad;
    double h=3.5*deg2rad;
    double xmin=-h/2.0,xmax=h/2.0;
    double ymin=0,ymax=w;
    double x,y,z;
    double pix=6.2/60.0/60.0*deg2rad;
    int numx,numy;
    numx=(xmax-xmin)/pix;
    numy=(ymax-ymin)/pix;
    cout<<numx<<"\t"<<numy<<endl;


    for(ii=0;ii<3;ii++) for(jj=0;jj<3;jj++) Amat[ii][jj]=rot[ii][jj];

    for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
    {
      x=xmin+pix*ii;
      y=ymin+pix*jj;
      z=0;

      Bmat[0][0]=x;
      Bmat[1][0]=y;
      Bmat[2][0]=z;

      AdotB(Amat,Bmat,Cmat,3,3,1);

      cout<<"("<<x<<","<<y<<","<<z<<")->("<<x2[0]+Cmat[0][0]<<","<<x2[1]+Cmat[1][0]<<","<<x2[2]+Cmat[2][0]<<")"<<endl;
    }

  }

  return 1;
}
