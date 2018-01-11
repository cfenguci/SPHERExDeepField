
#include "raj_code.h"
#include "global.h"
//double get_ZL(const int id_pts, const double x, const double y)
//{
//}

double ** bin_field(double **p)
{
  double ww=WIDTH_DETECTOR;
  double hh=HEIGHT_DETECTOR;
  double w=ww*deg2rad;
  double h=hh*deg2rad;

  double strip_width;
#ifdef _USE_OPENGL_
  strip_width=w/48.0*DISPLAY_REDUC;
#else
  strip_width=w/48.0;
#endif

  double xmin=-h/2.0,xmax=h/2.0;
  double ymin,ymax;
  get_strip(CURRENT_STRIP_ID,NUM_STRIP,ymin,ymax);

  double pix=dtr->pix;//big pixel
  int numx,numy;
  numx=(xmax-xmin)/pix;
  numy=(ymax-ymin)/pix;

  GRIDXY *field_coarse=new GRIDXY;
  field_coarse=NEW_GRIDXY(numx,numy,xmin,xmax,ymin,ymax);

  double **pLOW=create_2d_grid<double>(numx,numy);



  double pix0=pix/PIX_REDUC;
  int nx=(xmax-xmin)/pix0;
  int ny=(ymax-ymin)/pix0;
  cout<<"bin_field->high dims="<<nx<<"\t"<<ny<<endl;
  int i,j;
  int jj;
  for(i=0;i<nx;i++) for(j=0;j<ny;j++)
  {
    double x,y;
    x=xmin+pix0*i;
    y=ymin+pix0*j;
    if(j>=nx) jj=j-nx;
    else jj=j;
    double v=p[i][j];
    insert_grid(x,y,v,field_coarse);
  }

  flush_grid(field_coarse);

  for(i=0;i<numx;i++) for(j=0;j<numy;j++) pLOW[i][j]=field_coarse->grid[i][j];

  return pLOW;

}

void prepare_systematics()
{
  string filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/18424_IDrk.fits.txt";
  dt_dark_current=load_data_table(filename);
  cout<<"Syetematics->prepared "<<filename<<endl;
  //cout<<dt_dark_current->nrow<<"\t"<<dt_dark_current->ncol<<endl;



  //filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/thermalsig_1mK_eps.fits.txt";
  //filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/thermalsig_1mk_100s_eps_A.fits.txt";
  //filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/thermalsig_1mk_100s_eps_B.fits.txt";
  filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/OLD_thermalsig_1mK_eps.fits.txt";
  dt_thermo=load_data_table(filename);
  cout<<"Syetematics->prepared "<<filename<<endl;
  //cout<<dt_dark_current->nrow<<"\t"<<dt_dark_current->ncol<<endl;

  filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/zodi_A.txt";
  dt_zodi_A=load_data_table(filename);
  cout<<"Syetematics->prepared "<<filename<<endl;

  filename="/home/cfeng/computing/spherex/LVFscans/output/systematics/zodi_B.txt";
  dt_zodi_B=load_data_table(filename);
  cout<<"Syetematics->prepared "<<filename<<endl;


  int arraysize=2048;
  mat_dark_current=create_2d_grid<double>(arraysize,arraysize);
  mat_thermo=create_2d_grid<double>(arraysize,arraysize);
  mat_zodi_A=create_2d_grid<double>(arraysize,arraysize);
  mat_zodi_B=create_2d_grid<double>(arraysize,arraysize);


  int i,j;
  for(i=0;i<dt_dark_current->nrow;i++) 
  {
    int idx=dt_dark_current->table[i][0];
    int idy=dt_dark_current->table[i][1];
    mat_dark_current[idx][idy]=dt_dark_current->table[i][2];
    mat_thermo[idx][idy]=dt_thermo->table[i][2];
    mat_zodi_A[idx][idy]=dt_zodi_A->table[i][2];
    mat_zodi_B[idx][idy]=dt_zodi_B->table[i][2];
  }


  mat_dark_current_LOW=bin_field(mat_dark_current);
  mat_thermo_LOW=bin_field(mat_thermo);
  mat_zodi_A_LOW=bin_field(mat_zodi_A);
  mat_zodi_B_LOW=bin_field(mat_zodi_B);


}




double get_stamp_pixel(const double x, const double y, double **p)
{
  double ww=WIDTH_DETECTOR;
  double hh=HEIGHT_DETECTOR;
  double w=ww*deg2rad;
  double h=hh*deg2rad;

  double strip_width;
#ifdef _USE_OPENGL_
  strip_width=w/48.0*DISPLAY_REDUC;
#else
  strip_width=w/48.0;
#endif

  double xmin=-h/2.0,xmax=h/2.0;
  double ymin,ymax;
  get_strip(CURRENT_STRIP_ID,NUM_STRIP,ymin,ymax);

  double pix=dtr->pix/PIX_REDUC;//6.2''
  int numx,numy,arraysize=2048;
  numx=(xmax-xmin)/pix;
  numy=(ymax-ymin)/pix;

  int idx=(x-xmin)/pix;
  int idy=(y-ymin)/pix;
  //if(idy>=numx) idy-=numx;

  double v;

  if((idx>=0)&&(idx<numx)&&(idy>=0)&&(idy<numy)) v=p[idx][idy];
  else v=0;

  return v;
}


//double get_linear_thermo_ramp(x,y)
//{
//}

