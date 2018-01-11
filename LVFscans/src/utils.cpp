#include "lvf_scans.h"
/*
template<typename MATTYPE> MATTYPE **create_2d_grid(const int nx, const int ny)
{
  int i,j;
  MATTYPE **mat=(MATTYPE **)calloc(nx,sizeof(MATTYPE *));
  for(i=0;i<nx;i++) mat[i]=(MATTYPE *)calloc(ny,sizeof(MATTYPE));
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) mat[i][j]=(MATTYPE) MATTYPE_ZERO;
  return mat;
}

template<typename MATTYPE> void delete_2d_grid(MATTYPE ** mat, const int nx, const int ny)
{
  int i;
  for(i=0;i<nx;i++) free(mat[i]);
  free(mat);
}


template<typename MATTYPE> MATTYPE ***create_3d_grid(const int nx, const int ny, const int nz)
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
*/


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

void AdotB(double **A, double multi, const int dimA, const int dimB)
{
  int i,j;
  for(i=0;i<dimA;i++) for(j=0;j<dimB;j++) A[i][j]*=multi;

}



double det33(double **A)
{
return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-
A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+
A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

}


void rot_x(const double theta,double **xrot)
{
  xrot[0][0]=1;
  xrot[0][1]=0;
  xrot[0][2]=0;

  xrot[1][0]=0;
  xrot[1][1]=cos(theta);
  xrot[1][2]=-sin(theta);

  xrot[2][0]=0;
  xrot[2][1]=sin(theta);
  xrot[2][2]=cos(theta);

}


void rot_y(const double theta,double **yrot)
{
  yrot[0][0]=cos(theta);
  yrot[0][1]=0;
  yrot[0][2]=sin(theta);

  yrot[1][0]=0;
  yrot[1][1]=1;
  yrot[1][2]=0;

  yrot[2][0]=-sin(theta);
  yrot[2][1]=0;
  yrot[2][2]=cos(theta);

}


void rot_z(const double theta,double **zrot)
{
  zrot[0][0]=cos(theta);
  zrot[0][1]=-sin(theta);
  zrot[0][2]=0;

  zrot[1][0]=sin(theta);
  zrot[1][1]=cos(theta);
  zrot[1][2]=0;

  zrot[2][0]=0;
  zrot[2][1]=0;
  zrot[2][2]=1;

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

  int flag=0;
  if(modu==0) 
  {
    x[0]=0.0;;
    x[1]=0.0;
    x[2]=0.0;;

    flag=1;
    cout<<"Sigularity at: ("<<theta<<","<<phi<<")"<<endl;
  }
  else 
  {
    x[0]/=modu;
    x[1]/=modu;
    x[2]/=modu;

    flag=0;
  }
  return flag;
}


//for arbitrary rotation
void get_matW(double *u, double **matW)
{
  double ux,uy,uz;
  ux=u[0];
  uy=u[1];
  uz=u[2];
  
  matW[0][0]=0.0;matW[0][1]=-uz;matW[0][2]=uy;
  matW[1][0]=uz;matW[1][1]=0.0;matW[1][2]=-ux;
  matW[2][0]=-uy;matW[2][1]=ux;matW[2][2]=0.0;
}

void get_arbitrary_rot(const double phi, double *u, double **mat)
{
  
  int i,j;
  int dimA,dimB;
  dimA=3;
  dimB=3;
  double f1=sin(phi);
  double f2=2.0*sin(phi/2.0)*sin(phi/2.0);
  
  double **mat1=create_2d_grid<double>(3,3);
  double **mat2=create_2d_grid<double>(3,3);
  double **mat22=create_2d_grid<double>(3,3);
  double **mat3=create_2d_grid<double>(3,3);
  
  mat3[0][0]=1.0;mat3[0][1]=0.0;mat3[0][2]=0.0;
  mat3[1][0]=0.0;mat3[1][1]=1.0;mat3[1][2]=0.0;
  mat3[2][0]=0.0;mat3[2][1]=0.0;mat3[2][2]=1.0;
  
  
  get_matW(u,mat1);
  get_matW(u,mat2); AdotB(mat2,mat2,mat22,3,3,3);
  AdotB(mat1,f1,3,3);
  AdotB(mat22,f2,3,3);
  
  for(i=0;i<dimA;i++) for(j=0;j<dimB;j++)
  {
    mat3[i][j]+=mat1[i][j];
	mat3[i][j]+=mat22[i][j];
	mat[i][j]=mat3[i][j];
  }
  
  delete_2d_grid<double>(mat1,3,3);
  delete_2d_grid<double>(mat2,3,3);
  delete_2d_grid<double>(mat22,3,3);
  delete_2d_grid<double>(mat3,3,3);
}


double get_vec_len(double *vec)
{
  int i;
  double sum=0.0;
  for(i=0;i<3;i++) sum+=vec[i]*vec[i];
  return sqrt(sum);
}

void mat_cp(double **A,double **B)
{
  int i,j;
  for(i=0;i<3;i++) for(j=0;j<3;j++) B[i][j]=A[i][j];
}
void vec_cp(double *A,double *B)
{
for(int i=0;i<3;i++) B[i]=A[i];
}


COORD **create_2d_grid_coord(const int nx, const int ny)
{
  int i,j;
  COORD **mat=(COORD **)calloc(nx,sizeof(COORD *));
  for(i=0;i<nx;i++) mat[i]=(COORD *)calloc(ny,sizeof(COORD));
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) {mat[i][j].ra=0;mat[i][j].dec=0;}
  return mat;
}


void delete_2d_grid_coord(COORD ** mat, const int nx, const int ny)
{
  int i;
  for(i=0;i<nx;i++) free(mat[i]);
  free(mat);
}


COORD ***create_3d_grid_coord(const int nx, const int ny, const int nz)
{
  int i,j,k;
  COORD *** mat=(COORD ***)calloc(nx,sizeof(COORD **));
  for(i=0;i<nx;i++)
  {
    mat[i]=(COORD **)calloc(ny,sizeof(COORD *));
    for(j=0;j<ny;j++) mat[i][j]=(COORD *)calloc(nz,sizeof(COORD));
  }
  for(i=0;i<nx;i++)for(j=0;j<ny;j++)for(k=0;k<nz;k++) 
  {
    mat[i][j][k].ra=0.0;
    mat[i][j][k].dec=0.0;
  }
  return mat;
}


extern DETECTOR_HISTORY *DH;
void insert_detector(DETECTOR *p)
{
  int i,j;
  DH->num_dtr++;

  int id=DH->num_dtr;
  DH->dtr_arr[id]->numx=p->numx;
  DH->dtr_arr[id]->numy=p->numy;
  
  for(i=0;i<p->numx;i++) for(j=0;j<p->numy;j++) 
  {
    DH->dtr_arr[id]->fov_image[i][j]=p->fov_image[i][j];
	DH->dtr_arr[id]->fov_coord3d_x[i][j]=p->fov_coord3d_x[i][j];
	DH->dtr_arr[id]->fov_coord3d_y[i][j]=p->fov_coord3d_y[i][j];
	DH->dtr_arr[id]->fov_coord3d_z[i][j]=p->fov_coord3d_z[i][j];
	
  }
  
  cout<<"-->Inserted Detector View "<<id<<endl;

  
}
