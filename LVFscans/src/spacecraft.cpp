#include "lvf_scans.h"

#include "global.h"

void show_spacecraft_z(const double angle_THETA, const double angle_DELTA, double **rot_LVF, double *xpos_LVF)
{
  int i,j;
  int numx=dtr->numx;
  int numy=dtr->numy;
  double xmin=dtr->xmin;
  double ymin=dtr->ymin;
  double pix=dtr->pix;

  double len,ratio,r,xmax,ymax;
  double x1,x2,y1,y2;

  double center_LVF[3],center_LVF_3d[3];
  double axis_x[3],axis_y[3];
  double axis_x_3d[3],axis_y_3d[3]; 
  
  double xFOV[4],yFOV[4],zFOV[4];
  double xLVF[4],yLVF[4],zLVF[4];
  
  int id_pointing, id_step;  
  double strip_width=dtr->strip_width;  
  double phi_strip;
  


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
  double z_height;
  
  
  id_pointing=_ID_STEP_/NUM_STEPPING;
  id_step=_ID_STEP_%NUM_STEPPING;


//define the parameters
  len=(10*numx-1)*pix;

//cout<<"show_spacecraft_z->"<<len<<endl;
  
  x1=xmin;
  x2=xmin+(numx-1)*pix;
  y1=ymin;
  y2=ymin+(numy-1)*pix; 

  ratio=len+1;
  z_height=2*1.07;
  r=(z_height+1.07)/1.07;


  
  
 
  center_LVF[0]=(x1+x2)/2.0;
  center_LVF[1]=(y1+y2)/2.0;
  center_LVF[2]=0.0;  

  axis_x[0]=center_LVF[0]+len;
  axis_x[1]=center_LVF[1];
  axis_x[2]=0.0;
  
  axis_y[0]=center_LVF[0];
  axis_y[1]=center_LVF[1]+len;
  axis_y[2]=0; 
  


  xFOV[0]=x1;
  yFOV[0]=y1;
  zFOV[0]=0;

  xFOV[1]=x2;
  yFOV[1]=y1;
  zFOV[1]=0;

  xFOV[2]=x2;
  yFOV[2]=y2;
  zFOV[2]=0;

  xFOV[3]=x1;
  yFOV[3]=y2;
  zFOV[3]=0;  

  
  


 



  if(_ID_STEP_%NUM_STEPPING==0)
  //if(_ID_STEP_/NUM_STEPPING==0)
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

  AdotB(last_mat, mat_tmp1,mat_SKY_pointing,3,3,3);

  rot_z(id_step*strip_width*deg2rad,mat_slews);





//the camera
  get_normal_plane_axis(rot_LVF,xpos_LVF,
                        center_LVF,axis_x,axis_y,
                        center_LVF_3d,axis_x_3d,axis_y_3d);
#ifdef _reuse_
//for(i=0;i<3;i++) cout<<"Test->"<<xpos_LVF[i]<<"\t"<<center_LVF_3d[i]<<endl;

#ifdef _USE_OPENGL_
  //SHOW_normal_plane_axis(len,center_LVF_3d,axis_x_3d,axis_y_3d);  
#endif
  
  get_instrument(xFOV,yFOV,zFOV,xpos_LVF,rot_LVF);






//the detector
  get_fov_plane_axis(z_height, mat_LVF_pointing, xpos_LVF,
                     axis_spacecraft_x,axis_spacecraft_y,axis_spacecraft_z);  

#ifdef _USE_OPENGL_
  //SHOW_fov_plane_axis(xpos_LVF,axis_spacecraft_x,axis_spacecraft_y,axis_spacecraft_z); 
#endif

//get FOV center coordinates (theta, phi) for the first step
  get_FOV_center(axis_spacecraft_z, theta_camera_pointing, phi_camera_pointing);
  if(_ID_STEP_%NUM_STEPPING==0)
  {
    for(i=0;i<3;i++) last_camera_pointing[i]=axis_spacecraft_z[i][0];
    last_camera_theta=theta_camera_pointing;
    last_camera_phi=phi_camera_pointing;
  }
  
#ifdef _USE_OPENGL_
//show the slew step for the current wave
  //SHOW_vh_circles(strip_width,CURRENT_STRIP_ID,last_camera_phi,axis_spacecraft_z);
//show the slew steps for the first wave
  //SHOW_vh_circles(strip_width,0,last_camera_phi,axis_spacecraft_z);
#endif
  
#ifdef _USE_OPENGL_
  //SHOW_camera(r, xFOV, yFOV, zFOV, xpos_LVF, mat_SKY_pointing, mat_slews);
#endif
  
  record_detector_pixel(r,xpos_LVF,mat_SKY_pointing,mat_slews);
#endif//reuse

  pair_mat *mat_correction=instrument_maneuver();
  instrument(center_LVF_3d, mat_SKY_pointing,mat_correction);
  delete_maneuver_mat(mat_correction);
  
 
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


void get_normal_plane_axis(double **rot_LVF, double *xpos_LVF, double center_LVF[3],double axis_x[3], double axis_y[3],
double center_LVF_3d[3],double axis_x_3d[3], double axis_y_3d[3])
{
  int i,j;
  double **Amat=create_2d_grid<double>(3,3);
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);
  
  for(i=0;i<3;i++) for(j=0;j<3;j++) Amat[i][j]=rot_LVF[i][j];
  for(i=0;i<3;i++) Bmat[i][0]=center_LVF[i];

  AdotB(Amat,Bmat,Cmat,3,3,1);

  for(i=0;i<3;i++) center_LVF_3d[i]=xpos_LVF[i]+Cmat[i][0];

//in-plane x
  for(i=0;i<3;i++) Bmat[i][0]=axis_x[i];
  AdotB(Amat,Bmat,Cmat,3,3,1);
  for(i=0;i<3;i++) axis_x_3d[i]=xpos_LVF[i]+Cmat[i][0];

  for(i=0;i<3;i++) Bmat[i][0]=axis_y[i];
  AdotB(Amat,Bmat,Cmat,3,3,1);
  for(i=0;i<3;i++) axis_y_3d[i]=xpos_LVF[i]+Cmat[i][0];
  
  delete_2d_grid<double>(Amat,3,3);
  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);
}

void get_FOV_center(double **p, double &theta_camera_pointing, double &phi_camera_pointing)
{
  
  double radius=sqrt(p[0][0]*p[0][0]+p[1][0]*p[1][0]);
  double radius3d=sqrt(radius*radius+p[2][0]*p[2][0]);
  theta_camera_pointing=atan(radius/p[2][0]);
  theta_camera_pointing/=deg2rad;
  if(theta_camera_pointing<0) theta_camera_pointing+=180;

  phi_camera_pointing=atan2(p[1][0],p[0][0]);
  phi_camera_pointing/=deg2rad;
  if(phi_camera_pointing<0) phi_camera_pointing+=360.0;
}



void unit_vec(double p[3],double up[3])
{
  double mod=0.0;
  int i;
  mod=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  mod=sqrt(mod);
  if(mod<=0) 
  {
    cout<<"Error: unit_vec; mod zero"<<endl;
    throw 0;
  }
  for(i=0;i<3;i++) up[i]=p[i]/mod;
}

void SHOW_normal_plane_axis(const double len, double pz[3],double px[3], double py[3])
{
  int i;

  double upz[3],upx[3],upy[3];
  double pz1[3],px1[3],py1[3];
  double pz2[3],px2[3],py2[3];

  for(i=0;i<3;i++)
  {
    pz1[i]=pz[i];
    px1[i]=px[i]-pz[i];
    py1[i]=py[i]-pz[i];
  }
  for(i=0;i<3;i++)
  {
    unit_vec(pz1,upz);
    unit_vec(px1,upx);
    unit_vec(py1,upy);
  }
  for(i=0;i<3;i++)
  {
    pz1[i]=pz[i];
    px1[i]=px[i];
    py1[i]=py[i];
  }
  for(i=0;i<3;i++) 
  {
    pz2[i]=pz1[i]+upz[i]*len;
    px2[i]=px1[i]+upx[i]*len;
    py2[i]=py1[i]+upy[i]*len;
  }


//z
  glBegin(GL_LINES);
  glColor3f(0,0,1);
  glVertex3f(pz1[0],pz1[1],pz1[2]);
  glVertex3f(pz2[0],pz2[1],pz2[2]);
  glEnd();

//x
  glBegin(GL_LINES);
  glColor3f(1.0,0,0);
  glVertex3f(pz1[0],pz1[1],pz1[2]);
  glVertex3f(px2[0],px2[1],px2[2]);
  glEnd();
 
//y
  glBegin(GL_LINES);
  glColor3f(0,1,0);
  glVertex3f(pz1[0],pz1[1],pz1[2]);
  glVertex3f(py2[0],py2[1],py2[2]);
  glEnd();

}

void get_normal_plane_trans_mat(const double ratio, double center_LVF_3d[3],double axis_x_3d[3], double axis_y_3d[3],double **mat_normal)
{
  int i,j;
  double diff_x[3],diff_y[3],diff_z[3];
  double modx,mody,modz;
	
  for(i=0;i<3;i++)
  {
    diff_x[i]=axis_x_3d[i]-center_LVF_3d[i];
    diff_y[i]=axis_y_3d[i]-center_LVF_3d[i];
    diff_z[i]=(ratio-1.0)*center_LVF_3d[i];
  }
  modx=get_vec_len(diff_x);
  mody=get_vec_len(diff_y);
  modz=get_vec_len(diff_z);
  for(i=0;i<3;i++) 
  {
    diff_x[i]/=modx;
    diff_y[i]/=mody;
    diff_z[i]/=modz;
  }

  for(i=0;i<3;i++)
  {
    mat_normal[i][0]=diff_x[i];
    mat_normal[i][1]=diff_y[i];
    mat_normal[i][2]=diff_z[i];
  }
}

void get_fov_plane_axis(const double z_height, double **mat,double *pos,
                        double **px,double **py,double **pz)
{
  int i,j;
  double **Bmat=create_2d_grid<double>(3,1);
  double subcoordz[3],subcoordx[3],subcoordy[3];  

  int numx=dtr->numx;
  int numy=dtr->numy;
  double xmin=dtr->xmin;
  double ymin=dtr->ymin;
  double pix=dtr->pix;
  double len=(10*numx-1)*pix;
  
  subcoordx[0]=5*len;
  subcoordx[1]=0;
  subcoordx[2]=0;

  subcoordy[0]=0;
  subcoordy[1]=5*len;
  subcoordy[2]=0;  
  
  subcoordz[0]=0;
  subcoordz[1]=0;
  subcoordz[2]=z_height;
  
  for(i=0;i<3;i++) Bmat[i][0]=subcoordx[i];  
  AdotB(mat,Bmat,px,3,3,1);
  
  for(i=0;i<3;i++) Bmat[i][0]=subcoordy[i];
  AdotB(mat,Bmat,py,3,3,1);
  
  for(i=0;i<3;i++) Bmat[i][0]=subcoordz[i];
  AdotB(mat,Bmat,pz,3,3,1);
  
  for(i=0;i<3;i++)
  {
    px[i][0]+=pos[i];
	py[i][0]+=pos[i];
	pz[i][0]+=pos[i];
  }
  delete_2d_grid<double>(Bmat,3,1);
}

void SHOW_fov_plane_axis(double center_LVF_3d[3],double **axis_spacecraft_x,double **axis_spacecraft_y,double **axis_spacecraft_z)
{
  glLineWidth(3);
  glPushMatrix();
  glBegin(GL_LINES);
  glColor3f(0,0,1);
  glVertex3f(center_LVF_3d[0],center_LVF_3d[1],center_LVF_3d[2]);
  glVertex3f(axis_spacecraft_z[0][0],axis_spacecraft_z[1][0],axis_spacecraft_z[2][0]);
  glEnd();
  glPopMatrix();
  
  glPushMatrix();
  glBegin(GL_LINES);
  glColor3f(1,0,0);
  glVertex3f(center_LVF_3d[0],center_LVF_3d[1],center_LVF_3d[2]);
  glVertex3f(axis_spacecraft_x[0][0],axis_spacecraft_x[1][0],axis_spacecraft_x[2][0]);
  glEnd();
  glPopMatrix();
  
  glPushMatrix();
  glBegin(GL_LINES);
  glColor3f(0,1,0);
  glVertex3f(center_LVF_3d[0],center_LVF_3d[1],center_LVF_3d[2]);
  glVertex3f(axis_spacecraft_y[0][0],axis_spacecraft_y[1][0],axis_spacecraft_y[2][0]);
  glEnd();
  glPopMatrix();
  glLineWidth(1);
}

void SHOW_vh_circles(const double strip_width, const int current_strip, const double last_camera_phi, double **axis_spacecraft_z)
{
  int i;
  double radius=sqrt(axis_spacecraft_z[0][0]*axis_spacecraft_z[0][0]+axis_spacecraft_z[1][0]*axis_spacecraft_z[1][0]);
  double radius3d=sqrt(radius*radius+axis_spacecraft_z[2][0]*axis_spacecraft_z[2][0]);

  for(i=0;i<NUM_STEPPING;i++)
  {
    double phi_strip=last_camera_phi+current_strip*strip_width+i*strip_width;
    make_circle(radius3d,phi_strip,360,0.5,1.0,0);
  }  
  make_circle_h(radius3d,last_camera_theta,360,1,1,1);
}

void get_instrument(double xFOV[4],double yFOV[4],double zFOV[4], 
                  double *xpos_LVF, double **rot_LVF)
{
  int i,j,ii;
  double **Bmat=create_2d_grid<double>(3,1);
  double **BBmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);
  double **Xmat=create_2d_grid<double>(3,1);
	
  double camera_pointing[4][3];
  
  for(i=0;i<4;i++)
  {
    Bmat[0][0]=xFOV[i];
    Bmat[1][0]=yFOV[i];
    Bmat[2][0]=zFOV[i];

    AdotB(rot_LVF,Bmat,Cmat,3,3,1);

    for(ii=0;ii<3;ii++) camera_pointing[i][ii]=Cmat[ii][0]+xpos_LVF[ii];
  }  

/*
  double diffx[3],diffy[3],udiffx[3],udiffy[3],udiffz[3];
  for(i=0;i<3;i++)
  {
    diffx[i]=camera_pointing[0][i]-camera_pointing[1][i];
    diffy[i]=camera_pointing[0][i]-camera_pointing[3][i];
  }
  unit_vec(diffx,udiffx);
  unit_vec(diffy,udiffy);

  vecA_cross_vecB(udiffx,udiffy,udiffz);

  double hight_instrument[3];
  for(i=0;i<3;i++) hight_instrument[i]=hight*udiffz[i];
*/


  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(BBmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);
  delete_2d_grid<double>(Xmat,3,1);

  

#ifdef _USE_OPENGL_
  glPushMatrix();
  glColor3f(0.5,1,0.5);
  glBegin(GL_QUADS);
  for(i=0;i<4;i++) glVertex3f(camera_pointing[i][0],camera_pointing[i][1],camera_pointing[i][2]);
  glEnd();
  glPopMatrix();
#endif
  

}

void SHOW_camera(const double r, double xFOV[4],double yFOV[4],double zFOV[4], 
double *center_LVF, double **mat_SKY_pointing,double **mat_slews)
{
  int i,j,k,ii;
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);
  double **Xmat=create_2d_grid<double>(3,1);
	
  double camera_pointing_slews[4][3][1];
  double camera_pointing[4][3];
  
  for(i=0;i<4;i++)
  {
	Bmat[0][0]=r*xFOV[i];
    Bmat[1][0]=r*yFOV[i];
    Bmat[2][0]=r*zFOV[i];

	//AdotB(mat_slews,mat_spacecraft_z,mat_spacecraft_z_slews,3,3,3);
    AdotB(mat_SKY_pointing, Bmat,Cmat,3,3,1);

    for(ii=0;ii<3;ii++) 
	{
		camera_pointing[i][ii]=Cmat[ii][0]+last_camera_pointing[ii];
		Bmat[ii][0]=camera_pointing[i][ii];
	}
	AdotB(mat_slews,Bmat,Xmat,3,3,1);
	for(ii=0;ii<3;ii++) camera_pointing_slews[i][ii][0]=Xmat[ii][0];
  }  
  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);
  delete_2d_grid<double>(Xmat,3,1);

  

  


  double p1[3],p2[3],p3[3],p4[3];
  double alpha=0.4;
  reset_color();
  for(k=0;k<3;k++) 
  {
    p1[k]=camera_pointing_slews[0][k][0];
    p2[k]=camera_pointing_slews[1][k][0];
    p3[k]=camera_pointing_slews[2][k][0];
    p4[k]=camera_pointing_slews[3][k][0];
  }
  glPushMatrix();
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0,1.0,1.0,alpha);
  glBegin(GL_TRIANGLES);	
  glVertex3f(center_LVF[0],center_LVF[1],center_LVF[2]);
  glVertex3f(p1[0],p1[1],p1[2]);
  glVertex3f(p2[0],p2[1],p2[2]);
  glEnd();
  glBegin(GL_TRIANGLES);	
  glVertex3f(center_LVF[0],center_LVF[1],center_LVF[2]);
  glVertex3f(p2[0],p2[1],p2[2]);
  glVertex3f(p3[0],p3[1],p3[2]);
  glEnd();
  glBegin(GL_TRIANGLES);	
  glVertex3f(center_LVF[0],center_LVF[1],center_LVF[2]);
  glVertex3f(p3[0],p3[1],p3[2]);
  glVertex3f(p4[0],p4[1],p4[2]);
  glEnd();
  glBegin(GL_TRIANGLES);	
  glVertex3f(center_LVF[0],center_LVF[1],center_LVF[2]);
  glVertex3f(p4[0],p4[1],p4[2]);
  glVertex3f(p1[0],p1[1],p1[2]);
  glEnd();
  glPopMatrix();


}


void get_detector_pixel(const double r, double rFOV[3], 
double *xpos_LVF, double **mat_SKY_pointing,double **mat_slews, double rFOV_sky[3])
{
  int i,j,ii;
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);
  double **Xmat=create_2d_grid<double>(3,1);
	
  

  Bmat[0][0]=r*rFOV[0];
  Bmat[1][0]=r*rFOV[1];
  Bmat[2][0]=r*rFOV[2];

  AdotB(mat_SKY_pointing, Bmat,Cmat,3,3,1);

  for(ii=0;ii<3;ii++) Bmat[ii][0]=Cmat[ii][0]+last_camera_pointing[ii];
  
  AdotB(mat_slews,Bmat,Xmat,3,3,1);
  
  for(ii=0;ii<3;ii++) rFOV_sky[ii]=Xmat[ii][0];
  
  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);
  delete_2d_grid<double>(Xmat,3,1);

  

  

}


void record_detector_pixel(const double r, double *xpos_LVF, double **mat_SKY_pointing,double **mat_slews)
{
  int i,j;
  int ii,jj;

  double signal,noise;
  double noise_level=10.0;
  double eta1,eta2,eta3,eta4;
  
  eta1=0;//_rng_.rand_gauss();
  eta2=1.0;//_rng_.rand_gauss();


  
  double FOV_theta,FOV_phi;

  int numx,numy;
  double xmin,ymin,pix;
  double x,y,z;
  double rFOV[3],rFOV_sky[3],rFOV_sky_unit[3];
  double dis;
  hpint64 nside,idpix,id_parent;

  nside=_nside_;

  numx=dtr->numx;
  numy=dtr->numy;
  xmin=dtr->xmin;
  ymin=dtr->ymin;
  pix=dtr->pix;

  for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
  {
	eta3=_rng_.rand_gauss();

    x=xmin+pix*ii;
    y=ymin+pix*jj;

    z=0;

    rFOV[0]=x;
    rFOV[1]=y;
    rFOV[2]=z;
	

    get_detector_pixel(r, rFOV,xpos_LVF, mat_SKY_pointing,mat_slews,rFOV_sky);
	
    dis=get_vec_len(rFOV_sky);


//cout<<"record_detector_pixel->"<<dis<<endl;

	if(dis!=0) for(i=0;i<3;i++) rFOV_sky_unit[i]=rFOV_sky[i]/dis;
	else for(i=0;i<3;i++) rFOV_sky_unit[i]=1e-8;
    spherical_2_cartesian(rFOV_sky_unit[0],rFOV_sky_unit[1],rFOV_sky_unit[2],FOV_theta,FOV_phi);

	//cout<<FOV_theta*deg2rad<<"\t"<<FOV_phi*deg2rad<<endl;
    ang2pix_ring64(nside,FOV_theta*deg2rad,FOV_phi*deg2rad,&idpix);

//idpix=1;

    id_parent=idpix;//get_index_mapping_kid2parent(idpix);

    signal=_map_GAL_[id_parent];
    noise=eta3*noise_level;

    _map_NAIVE_[id_parent]+=(signal+noise);
    _map_LVF_scan_[id_parent]++;

	double pix=PIX_REDUC*6.2/60.0/60.0*deg2rad;
/*
    dtr->fov_coord[ii][jj].ra=x+2*pix*eta1;
    dtr->fov_coord[ii][jj].dec=y+2*pix*eta1;
    dtr->fov_image[ii][jj]=eta2*(signal+noise);

    dtr->fov_coord3d_x[ii][jj]=rFOV_sky[0];
    dtr->fov_coord3d_y[ii][jj]=rFOV_sky[1];
    dtr->fov_coord3d_z[ii][jj]=rFOV_sky[2];
*/

  }

  double rr,gg,bb;
  double v;
  
  /*
#ifdef _USE_OPENGL_
  glPushMatrix();
  glBegin(GL_POINTS);
  for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
  {
    v=dtr->fov_image[ii][jj]*1e7;

    x=dtr->fov_coord3d_x[ii][jj];
    y=dtr->fov_coord3d_y[ii][jj];
    z=dtr->fov_coord3d_z[ii][jj];

    color2rgb(v,rr,gg,bb);
    glColor3f(rr,gg,bb);
    glVertex3f(x,y,z);

  }
  glEnd();
  glPopMatrix();
#endif
*/

#ifdef _SHOW_DECTECTOR_HISTORY_
  insert_detector(dtr);
  
  //cout<<DH->num_dtr<<endl;
  
  if(DH->num_dtr%30==0)
  {
  
        glPushMatrix();
      glBegin(GL_POINTS);
  
    for(int id_dtr=0;id_dtr<DH->num_dtr;id_dtr++)
    {
		//cout<<id_dtr<<endl;

      for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
      {
        v=DH->dtr_arr[id_dtr]->fov_image[ii][jj]*1e7;

        x=DH->dtr_arr[id_dtr]->fov_coord3d_x[ii][jj];
        y=DH->dtr_arr[id_dtr]->fov_coord3d_y[ii][jj];
        z=DH->dtr_arr[id_dtr]->fov_coord3d_z[ii][jj];

        color2rgb(v,rr,gg,bb);
        glColor3f(rr,gg,bb);
        glVertex3f(x,y,z);

      }

    }  
  	  glEnd();  
	  glPopMatrix();
	  
	  
	    make_movie(DH->num_dtr/30,"scan/patches",windW,windH);
	  
  }
  

  
  
#endif

}


pair_mat * instrument_maneuver()
{
  pair_mat *pm=new pair_mat;
  pm->mat_theta=create_2d_grid<double>(3,3);
  pm->mat_phi=create_2d_grid<double>(3,3);

  int id_pointing, id_step;
  double strip_width=dtr->strip_width;
  double dt=dtr->delta_t;

  id_pointing=_ID_STEP_/NUM_STEPPING;
  id_step=_ID_STEP_%NUM_STEPPING;

  double omega_theta=dtr->omega_theta;
  double delta_theta=0;//-id_step*omega_theta*dt;
  double delta_phi=id_step*strip_width/deg2rad;
#ifdef _USE_OPENGL_
  cout<<"instrument_maneuver->"<<strip_width/deg2rad<<endl;
  //cout<<"instrument_maneuver->"<<delta_theta<<"\t"<<delta_phi<<endl;
#endif

  rot_y(delta_theta*deg2rad,pm->mat_theta);
  rot_z(delta_phi*deg2rad,pm->mat_phi);

  return pm;
}

void delete_maneuver_mat(pair_mat *pm)
{
  delete_2d_grid<double>(pm->mat_theta,3,3);
  delete_2d_grid<double>(pm->mat_phi,3,3);
  delete pm;
}





