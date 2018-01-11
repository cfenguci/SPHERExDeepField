#include "lvf_scans.h"

extern DETECTOR *dtr;
extern int _ID_STEP_;

extern Healpix_Map<double> _map_GAL_;
extern Healpix_Map<double> _map_NAIVE_;
extern Healpix_Map<double> _map_LVF_scan_;

extern planck_rng _rng_;

extern hpint64 _nside_;

extern double **last_mat,*last_pos,*last_camera_pointing,last_camera_theta,last_camera_phi;

extern Healpix_Map<hpint64> _index_mapping_;


extern DETECTOR_HISTORY *DH;



void show_spacecraft_z_patch(const double angle_THETA, const double angle_DELTA, double **rot_LVF, double *xpos_LVF)
{
  int i,j;
  int numx=dtr->numx;
  int numy=dtr->numy;
  double xmin=dtr->xmin;
  double ymin=dtr->ymin;
  double pix=dtr->pix;
//show the spacecraft +z axis
  double len=(3*numx-1)*pix;
  
  double x1=xmin;
  double x2=xmin+(numx-1)*pix;
  double y1=ymin;
  double y2=ymin+(numy-1)*pix; 



  double center_LVF[3],center_LVF_3d[3];
  double axis_x[3],axis_y[3];
  double axis_x_3d[3],axis_y_3d[3]; 
  

  
  double xFOV[4],yFOV[4],zFOV[4];
  double xLVF[4],yLVF[4],zLVF[4];
  
  double ratio=len+1;
  
  int id_pointing, id_step;
  
  double w=WIDTH_DETECTOR;
  
  double strip_width=w/48*DISPLAY_REDUC;
  

  
  double phi_strip;
  
  double r=5;

  double **Amat=create_2d_grid<double>(3,3);
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);
  double **mat_tmp1=create_2d_grid<double>(3,3);
  
  double **mat_x=create_2d_grid<double>(3,3);
  double **mat_y=create_2d_grid<double>(3,3);
  
  double **mat_normal=create_2d_grid<double>(3,3);
  double **mat_LVF_pointing=create_2d_grid<double>(3,3);
  double **mat_SKY_pointing=create_2d_grid<double>(3,3);
  
  double **mat_slews=create_2d_grid<double>(3,3);
  double **mat_spacecraft_z_slews=create_2d_grid<double>(3,3);

  double **axis_spacecraft_x=create_2d_grid<double>(3,1);  
  double **axis_spacecraft_y=create_2d_grid<double>(3,1);
  double **axis_spacecraft_z=create_2d_grid<double>(3,1);
  
  double **fov_camera=create_2d_grid<double>(3,1);
  double theta_camera_pointing, phi_camera_pointing;
  
  
  id_pointing=_ID_STEP_/NUM_STEPPING;
  id_step=_ID_STEP_%NUM_STEPPING;
  
  
 
  center_LVF[0]=(x1+x2)/2.0;
  center_LVF[1]=(y1+y2)/2.0;
  center_LVF[2]=0.0;  

  axis_x[0]=center_LVF[0]+len;
  axis_x[1]=center_LVF[1];
  axis_x[2]=0.0;
  
  axis_y[0]=center_LVF[0];
  axis_y[1]=center_LVF[1]+len;
  axis_y[2]=0;
  
  
  
  


  xFOV[0]=xmin;
  yFOV[0]=ymin;
  zFOV[0]=0;

  xFOV[1]=xmin+(numx-1)*pix;
  yFOV[1]=ymin;
  zFOV[1]=0;

  xFOV[2]=xmin+(numx-1)*pix;
  yFOV[2]=ymin+(numy-1)*pix;
  zFOV[2]=0;

  xFOV[3]=xmin;
  yFOV[3]=ymin+(numy-1)*pix;
  zFOV[3]=0;  

  
  


 



  if(_ID_STEP_%NUM_STEPPING==0)
  {
    mat_cp(rot_LVF, last_mat);
    vec_cp(xpos_LVF, last_pos);
  }

//delta
  rot_y(angle_DELTA*deg2rad,mat_y);
//theta
  rot_x(angle_THETA*deg2rad,mat_x);

  AdotB(mat_y,mat_x,mat_tmp1,3,3,3);
  AdotB(rot_LVF,mat_tmp1,mat_LVF_pointing,3,3,3);
  //mat_cp(rot_LVF, mat_LVF_pointing);
  AdotB(last_mat, mat_tmp1,mat_SKY_pointing,3,3,3);

  rot_z(id_step*strip_width*deg2rad,mat_slews);

  get_normal_plane_axis(rot_LVF,xpos_LVF,center_LVF,axis_x,axis_y,center_LVF_3d,axis_x_3d,axis_y_3d);


  
  get_instrument(xFOV,yFOV,zFOV,xpos_LVF,rot_LVF);
  //get_normal_plane_trans_mat(ratio,center_LVF_3d,axis_x_3d,axis_y_3d,mat_normal);    

  get_fov_plane_axis(len, mat_LVF_pointing,xpos_LVF,axis_spacecraft_x,axis_spacecraft_y,axis_spacecraft_z);  



//get FOV center coordinates (theta, phi) for the first step
  get_FOV_center(axis_spacecraft_z, theta_camera_pointing, phi_camera_pointing);
  if(_ID_STEP_%NUM_STEPPING==0)
  {
    for(i=0;i<3;i++) last_camera_pointing[i]=axis_spacecraft_z[i][0];
    last_camera_theta=theta_camera_pointing;
    last_camera_phi=phi_camera_pointing;
  }
  

  

  
  record_detector_pixel(r,xpos_LVF,mat_SKY_pointing,mat_slews);
  
 
  delete_2d_grid<double>(Amat,3,3);
  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);
  delete_2d_grid<double>(mat_tmp1,3,3);
  
  delete_2d_grid<double>(mat_x,3,3);
  delete_2d_grid<double>(mat_y,3,3);
  
  delete_2d_grid<double>(mat_normal,3,3);
  delete_2d_grid<double>(mat_SKY_pointing,3,3);
  delete_2d_grid<double>(mat_LVF_pointing,3,3);
  
  delete_2d_grid<double>(mat_slews,3,3);
  delete_2d_grid<double>(mat_spacecraft_z_slews,3,3);
  
  delete_2d_grid<double>(axis_spacecraft_x,3,1);
  delete_2d_grid<double>(axis_spacecraft_y,3,1);  
  delete_2d_grid<double>(axis_spacecraft_z,3,1);

  delete_2d_grid<double>(fov_camera,3,1);
}



