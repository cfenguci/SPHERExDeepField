#include "lvf_scans.h"

#include "global.h"

void make_circle(const double radius, const double phi, const int num_grid, const double r, const double g, const double b)
{
  int i;
  double theta;
  double dtheta=360.0/num_grid;
  double x[3];
  //glPushMatrix();
  glBegin(GL_LINES);
  for(i=0;i<num_grid;i++)
  {
    theta=i*dtheta;
    x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
    x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
    x[2]=radius*cos(deg2rad*theta);
    //glColor3f(0.8,0.3,0.5);
    glColor3f(r,g,b);
    glVertex3f(x[0],x[1],x[2]);
  }
  glEnd();
  //glPopMatrix();
}

void make_circle_h(const double radius, const double theta, const int num_grid,const double r, const double g, const double b)
{
  int i;
  double phi;
  double dphi=360.0/num_grid;
  double x[3];
  //glPushMatrix();
  glBegin(GL_LINES);
  for(i=0;i<num_grid;i++)
  {
    phi=i*dphi;
    x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
    x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
    x[2]=radius*cos(deg2rad*theta);
    glColor3f(r,g,b);
    glVertex3f(x[0],x[1],x[2]);
  }
  glEnd();
  
  //glPopMatrix();
  //cout<<"radius="<<radius<<endl;
}

void make_circle_h(const double radius, const double theta, const int num_grid,COORD3D *ring)
{
  int i;
  double phi;
  double dphi=360.0/num_grid;
  double x[3];

  for(i=0;i<num_grid;i++)
  {
    phi=i*dphi;
    x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
    x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
    x[2]=radius*cos(deg2rad*theta);

    for(int k=0;k<3;k++) ring[i].r[k]=x[k];

  }

}

void make_circle_v(const double radius, const double phi, const int num_grid,COORD3D *ring)
{
  int i;
  double theta;
  double dtheta=360.0/num_grid;
  double x[3];

  for(i=0;i<num_grid;i++)
  {
    theta=i*dtheta;
    x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
    x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
    x[2]=radius*cos(deg2rad*theta);

    for(int k=0;k<3;k++) ring[i].r[k]=x[k];

  }

}



void generate_pointing(const double radius, const double theta_0, const double phi_0,const double omega_theta, const double omega_phi,const double t,
double *x,double *xn, 
double &theta_out, double &phi_out)
{
  double theta,phi;
  
  theta=theta_0+omega_theta*t;
  phi=phi_0+omega_phi*t;

  sphere_normal_vec(theta,phi,xn); 

  x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
  x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
  x[2]=radius*cos(deg2rad*theta);

  theta_out=theta;
  phi_out=phi;
} 

void generate_pointing_slews(const double radius, const double theta_0, const double phi_0, const double dt, const double omega_theta, 
const double omega_phi,const int id_slews, const double strip_width, double *x,double &theta_out, double &phi_out)
{
  double theta,phi;
  theta=theta_0+omega_theta*dt*id_slews;
  phi=phi_0+omega_phi*dt*id_slews+strip_width*id_slews;

  x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
  x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
  x[2]=radius*cos(deg2rad*theta);

  theta_out=theta;
  phi_out=phi;
}


//deg
void spherical_2_cartesian(const double x, const double y, const double z, double &theta, double &phi)
{
  double ratio1,ratio2;
  double r=sqrt(x*x+y*y+z*z);

  ratio1=z/r;

  theta=acos(ratio1)/deg2rad;

  phi=atan2(y,x)/deg2rad;
  if(y<0) phi+=360.0;

}

void reset_color()
{
  GLfloat white[3];
  white[0]=1;
  white[1]=1;
  white[2]=1;
  glColor3fv(white);
}

void draw_sphere(GLUquadric *tmp_quad)
{

  reset_color();

  double rot=_ID_STEP_*185.0/3600.0*360.0/24.0;

  glRotatef(rot, 0.0f, 0.0f, 1.0f);
  //glRotatef(-90,1.0f,0.0f,0.0f);
  

  gluSphere(tmp_quad,RADIUS_EARTH,100,100);

}

void earth_texture()
{
  draw_sphere(quad_earth);
  id_sphere++;
}


void draw_sky(const double r, const int num_phi, const int num_theta, COORD3D **sphere)
{


  //glPushMatrix();
  //glRotatef(-90,1.0f,0.0f,0.0f);

  int i,j;
  //int num_theta=100;
  //int num_phi=200;
  double dtheta=180.0/num_theta;
  double dphi=360.0/num_phi;
  double theta,phi;
  double x,y,z;

  for(i=0;i<num_phi;i++) for(j=0;j<num_theta;j++)
  { 
    //glBegin(GL_POINTS);
    theta=j*dtheta;
    phi=dphi*i;

    theta*=deg2rad;
    phi*=deg2rad;

    x=r*sin(theta)*cos(phi);
    y=r*sin(theta)*sin(phi);
    z=r*cos(theta);

    sphere[i][j].r[0]=x;
    sphere[i][j].r[1]=y;
    sphere[i][j].r[2]=z;

    //cout<<"Sphere->"<<i<<"\t"<<j<<endl;

 
    //glColor3f(0.1,0.1,0.2);
    //glVertex3f(x,y,z);     
    
    //glEnd();
  }

  //glutWireSphere(3.0, 10,20);
  //glPopMatrix();



}






void get_LVF_image(const int idi,const int idj,const int idpix, const double dec, const double ra, DETECTOR *dt)
{
  dt->fov_image[idi][idj]=_map_GAL_[idpix];
  dt->fov_coord[idi][idj].dec=dec;
  dt->fov_coord[idi][idj].ra=ra;

  if(_max_ra_<=ra) _max_ra_=ra;
  if(_max_dec_<=dec) _max_dec_=dec;
  if(_min_ra_>=ra) _min_ra_=ra;
  if(_min_dec_>=dec) _min_dec_=dec;


}




void show_sub_window()
{
  int i;
  double x[4],y[4],z[4];

 

  int ii,jj;


  glBegin(GL_POINTS);
  for(jj=0;jj<dtr->numy;jj++)
  for(ii=0;ii<dtr->numx;ii++)
  {
    double rr,gg,bb;
    color2rgb(dtr->fov_image[ii][jj]*1e7,rr,gg,bb);

    glColor3f(rr,gg,bb);

    double x=dtr->fov_coord[ii][jj].dec/deg2rad-WIDTH_DETECTOR/4.0;
    double y=dtr->fov_coord[ii][jj].ra/deg2rad;


    glVertex3f(x,y,0);
  }
  glEnd();

  cout<<"Sub window->"<<_ID_STEP_<<endl;
//glPopMatrix();
}


void FOV_image(const int RecordSky, DETECTOR *dt,double **rot_LVF, double *xpos_LVF)
{
  show_spacecraft_z(TEST_THETA,TEST_NOD,rot_LVF,xpos_LVF);
  //show_spacecraft_z_patch(5.0,35.0,rot_LVF,xpos_LVF);
}



void init()
{
  for(int i=0;i<_map_LVF_scan_.Npix();i++) _map_LVF_scan_[i]=0.0;
  for(int i=0;i<_map_NAIVE_.Npix();i++) _map_NAIVE_[i]=0.0;
  
  CAMERA_X=DefaultCameraX;
  CAMERA_Y=DefaultCameraY;
  CAMERA_Z=DefaultCameraZ;    

  rotate_background();

  LVF_track=new_LVF_track(2*NUM_STEPPING,360);

  Count_NCP_pointings=0;
  Count_NCP_frames=0;

  INIT_NCP_flat();
}


void get_position_sub(void)
{

  glClear(GL_COLOR_BUFFER_BIT );
  glLoadIdentity();
  gluLookAt(0,0,-9,0,0,0,0,1,0);
  show_sub_window();

  glutSwapBuffers();

  glutPostRedisplay();


  make_movie("fov/frame",200,400);

}

void update_sat_theta()
{
_ID_STEP_=66;
TEST_THETA-=1;
}

void update_sat_delta()
{
_ID_STEP_=66;
TEST_NOD+=1;
}


void update_cant_angle()
{
  double delta_phi=_GLOBAL_CLOCK_;
  //if(TEST_CANT>0) TEST_CANT-=delta_phi;
  //if(TEST_CANT<=0) TEST_CANT+=delta_phi;
//  TEST_CANT-=delta_phi;

TEST_CANT-=1.0;
//if(TEST_CANT<10) TEST_CANT+=1.0;
//if(TEST_CANT>10) TEST_CANT-=1.0;
//if(TEST_CANT>70) TEST_CANT-=1.0;
  _ID_STEP_=66;
}


void update_CAMERA()
{
  //_ID_STEP_=11;

  double BOUND=1;
  double multiply=1.0+0.1*sqrt(_GLOBAL_CLOCK_);
  int step=0;
  
  /*
  int flag_zoom_in=0;
  if(multiply>=BOUND) 
  {
    multiply=BOUND;
    step=1;
	_FLAG_show_axis_=0;
	
	flag_zoom_in=1;
  }
  */

  double default_x=DefaultCameraX,default_y=DefaultCameraY,default_z=DefaultCameraZ;
  
//  CAMERA_X=default_x*multiply;
//  CAMERA_Y=default_y*multiply;
//  CAMERA_Z=default_z*multiply;


  
  if((_GLOBAL_CLOCK_>0)&&(_GLOBAL_CLOCK_<30))
  {
    CAMERA_X=default_x*multiply;
    CAMERA_Y=default_y*multiply;
    CAMERA_Z=default_z*multiply;
    step=0;
    double r=BOUND*sqrt(default_x*default_x+default_z*default_z);
    double rot=(_GLOBAL_CLOCK_)*2.0*deg2rad;	  

    CAMERA_X=r*cos(rot);
    CAMERA_Z=r*sin(rot);
	  //cout<<rot/deg2rad<<"\t"<<CAMERA_X<<"\t"<<CAMERA_Z<<endl;
  }
  
  if((_GLOBAL_CLOCK_>=30)&&(_GLOBAL_CLOCK_<43))
  {
      //multiply=1.0+_GLOBAL_CLOCK_;
    CAMERA_Z+=(0.1*(_GLOBAL_CLOCK_-30));
  }
  if(_GLOBAL_CLOCK_>=43)
  {
    FLAG_show_ncp=0;
    _ID_STEP_++;
    _CLOCK_SUN_++;
  }
  

	
  cout<<"Get the clock sequence -> "<<_GLOBAL_CLOCK_<<"\t"<<_ID_STEP_<<"\t"<<multiply<<endl;
}



void get_strip(const int id_strip, const int num, double &ymin, double &ymax)
{

  double strip_width;
  double w=WIDTH_DETECTOR;
#ifdef _USE_OPENGL_
  strip_width=w/48.0*DISPLAY_REDUC;
#else
  strip_width=w/48.0;
#endif

  strip_width*=deg2rad;

  double half=w/2.0*deg2rad;
  //ymin=-WIDTH_DETECTOR/2.0*deg2rad+id_strip*strip_width;
  ymin=id_strip*strip_width-half;
  ymax=ymin+num*strip_width;
}



DETECTOR * init_detector_view(const double ww,const double hh)
{
  double w=ww*deg2rad;
  double h=hh*deg2rad;
  
  double strip_width;
#ifdef _USE_OPENGL_
  strip_width=w/48.0*DISPLAY_REDUC;
#else
  strip_width=w/48.0;
#endif

  double xmin=-h/2.0,xmax=h/2.0;
  //double ymin=0,ymax=strip_width;
  double ymin,ymax;
  get_strip(CURRENT_STRIP_ID,NUM_STRIP,ymin,ymax);

  double pix=PIX_REDUC*6.15234375/60.0/60.0*deg2rad;
  int numx,numy;
  numx=(xmax-xmin)/pix;
  numy=(ymax-ymin)/pix;
  
  //cout<<"Detector dimensions: "<<numx<<"\t"<<numy<<endl;


  DETECTOR *dt=new DETECTOR;
  dt->numx=numx;
  dt->numy=numy;
  dt->xmin=xmin;
  dt->ymin=ymin;
  dt->pix=pix;
  dt->fov_image=create_2d_grid<double>(dt->numx+1,dt->numy+1);

  dt->fov_coord=create_2d_grid_coord(dt->numx+1,dt->numy+1);

  dt->fov_coord3d_x=create_2d_grid<double>(dt->numx+1,dt->numy+1);
  dt->fov_coord3d_y=create_2d_grid<double>(dt->numx+1,dt->numy+1);
  dt->fov_coord3d_z=create_2d_grid<double>(dt->numx+1,dt->numy+1);

  dt->strip_width=strip_width;


  //cout<<"Detector view prepared!"<<endl;
  
  return dt;
}

void delete_detector_view()
{
  delete_2d_grid<double>(dtr->fov_image,dtr->numx+1,dtr->numy+1);
  delete_2d_grid_coord(dtr->fov_coord,dtr->numx+1,dtr->numy+1);
  delete_2d_grid<double>(dtr->fov_coord3d_x,dtr->numx+1,dtr->numy+1);
  delete_2d_grid<double>(dtr->fov_coord3d_y,dtr->numx+1,dtr->numy+1);
  delete_2d_grid<double>(dtr->fov_coord3d_z,dtr->numx+1,dtr->numy+1);
  delete dtr;
}




void get_pole(const int id)
{
  int i,j;
  string file;
  //file="/home/cmb/computing/spherex/VLFscans/output/VLF_scan.fits";
  hpint64 nside=_map_LVF_scan_.Nside(),idpix;
  //Healpix_Map<double> map(nside,RING,SET_NSIDE);
  //read_Healpix_map_from_fits(file,map,1,2);
  vec3 vec;

  double xmin=-1.0,xmax=1.0;
  double ymin=-1.0,ymax=1.0;
  double dx=0.01,dy=0.01;
  int numx=(xmax-xmin)/dx;
  int numy=(ymax-ymin)/dy;
  double x,y;

  double **pole_count=create_2d_grid<double>(numx,numy);
  double **pole_count0=create_2d_grid<double>(numx,numy);

  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    pole_count[i][j]=0.0;
    pole_count0[i][j]=0.0;
  }

  double theta,phi;
  for(i=0;i<_map_LVF_scan_.Npix();i++)
  {
    pix2ang_ring64(nside, i, &theta, &phi);
    theta/=deg2rad;
    phi/=deg2rad;

    if(theta<=45.0)
    {
      vec=_map_LVF_scan_.pix2vec(i);
      int idx=(vec.x-xmin)/dx;
      int idy=(vec.y-ymin)/dy;
      pole_count[idx][idy]+=_map_LVF_scan_[i];
      pole_count0[idx][idy]++;
    }

  }


  ofstream out;
  string null="";
  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%d",id);
  cout<<id<<endl;
  out.open((null+"/home/cfeng/computing/spherex/LVFscans/output/pole_map_"+seg+".txt").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    if(pole_count0[i][j]<=0) out<<i<<"\t"<<j<<"\t"<<0<<endl;
    else out<<i<<"\t"<<j<<"\t"<<pole_count[i][j]<<endl;
  }
  out.close();


}

void timer(int v) 
{
  glutPostRedisplay();
  glutTimerFunc(100, timer, v);
}


void timer1(int v)
{
  glutPostRedisplay();
  glutTimerFunc(100, timer1, v);
}


void reshape(GLint w, GLint h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, GLfloat(w) / GLfloat(h), 1.0, -5.0);
  glMatrixMode(GL_MODELVIEW);
}


void resize_sub(int width, int height)
{
    const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}



#define ESC 27
float row, xpos = -10.0 , ypos = -10.0 , zpos = 18.0;// initial position
float deltaMove = 0.0; // initially camera doesn't move

//Camera direction
float lx = 1.0, ly = 1.0, lz = 0.0; //initially points along y-axis
float angle = 0.0; // angle of rotation for the camera direction
float anglez = 0.0;
float deltaAngle = 0.0; //additional angle change when dragging
float deltaAngley = 0.0;

//Mouse drag control
int isDragging = 0; //true when dragging
int xDragStart = 0; //records the x-coordinate when dragging starts
int yDragStart = 0; //records the y-coordinate when dragging starts



void mouseMove(int xpos, int ypos)
{
  if (isDragging) 
  { 		
    deltaAngle = (xpos - xDragStart) * 0.005;
    deltaAngley = (ypos - yDragStart) * 0.005;
    lx = sin(angle + deltaAngle);
    ly = -cos(angle + deltaAngle);
    lz = -sin(anglez + deltaAngley);
  }
}

void mouseButton(int button, int state, int xpos, int ypos) 
{
  if (button == GLUT_LEFT_BUTTON) 
  { 
    if (state == GLUT_DOWN) 
    { 			
      isDragging = 1; 
      xDragStart = xpos; 
      yDragStart = ypos;
      id_sphere++;
//cout<<xpos<<"\t"<<ypos<<endl;
      glutPostRedisplay();
    } 
    else  
    { 
      angle += deltaAngle; 
      anglez += deltaAngley;
      isDragging = 0; 
    }
  }
}

void make_movie(const int id, string name_movie, const int w, const int h)
{
  GLint components;
  GLint width, height;
  components=3;
  width=w;
  height=h;

  GLubyte *data = (GLubyte *)malloc(components * width * height*sizeof(GLubyte));
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);


  int i,j;
  char seg[1024];
  ofstream out;

  sprintf(seg,"%04d",id);
  char fname[1024];

  memset(fname,0,1024);
  strcat(fname,"/home/cfeng/computing/spherex/LVFscans/output/framedata/");
  strcat(fname,name_movie.c_str());
  strcat(fname,"_");
  strcat(fname,seg);
  strcat(fname,".bmp");
  cout<<fname<<endl;
  int cnt=0;
  WriteBMP(fname, width, height, data);
  free(data);
}

void make_movie(string name_movie, const int w, const int h)
{
  GLint components;
  GLint width, height;
  components=3;
  width=w;
  height=h;

  GLubyte *data = (GLubyte *)malloc(components * width * height*sizeof(GLubyte));
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);


  int i,j;
  char seg[1024];
  ofstream out;

  sprintf(seg,"%04d",_GLOBAL_CLOCK_);
  char fname[1024];

  memset(fname,0,1024);
  strcat(fname,"/home/cfeng/computing/spherex/LVFscans/output/framedata/");
  strcat(fname,name_movie.c_str());
  strcat(fname,"_");
  strcat(fname,seg);
  strcat(fname,".bmp");
  cout<<fname<<endl;
  int cnt=0;
  WriteBMP(fname, width, height, data);
  free(data);
}


void load_fits_map_cxx(string file)
{
  read_Healpix_map_from_fits(file,_map_GAL_,1,2);
  tsize nmod=_map_GAL_.replaceUndefWith0();
  if (nmod!=0) cout << "WARNING: replaced " << nmod << " undefined map pixels with a value of 0" << endl;
  if (_map_GAL_.Scheme()==NEST) _map_GAL_.swap_scheme();
  cout<<"Loaded "<<file<<endl;
}

void load_fits_map_cxx(string file, Healpix_Map<double> &HealpixDATA)
{
  read_Healpix_map_from_fits(file,HealpixDATA,1,2);
  tsize nmod=HealpixDATA.replaceUndefWith0();
  if (nmod!=0) cout << "WARNING: replaced " << nmod << " undefined map pixels with a value of 0" << endl;
  if (HealpixDATA.Scheme()==NEST) HealpixDATA.swap_scheme();
  cout<<"Loaded "<<file<<endl;
}

void write_grid_data(DETECTOR *dtr)
{
  ofstream out;
  int i,j;
  string null="";
  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%d",_ID_FOV_record_);
  cout<<"_ID_FOV_record_"<<"\t"<<seg<<endl;

  out.open((null+"/home/cfeng/computing/spherex/LVFscans/output/fov/FOV_Image_"+seg+".txt").c_str());
  for(i=0;i<dtr->numx;i++) for(j=0;j<dtr->numy;j++) out<<i<<"\t"<<j<<"\t"<<dtr->fov_coord[i][j].ra<<"\t"<<dtr->fov_coord[i][j].dec<<"\t"<<dtr->fov_image[i][j]<<endl;
  out.close();


  //cout<<_min_ra_<<"\t"<<_max_ra_<<endl;
  //cout<<_min_dec_<<"\t"<<_max_dec_<<endl;

  _ID_FOV_record_++;


}






void postprocessing()
{
  string outfile;
  int i;
  for(i=0;i<_map_NAIVE_.Npix();i++) 
  {
    double nhits=_map_LVF_scan_[i];
    if(nhits!=0) _map_NAIVE_[i]/=nhits;
    else _map_NAIVE_[i]=0.0;
  }


  get_naive_map_powspec();


  outfile="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/naive_scan.fits";
  write_Healpix_map_to_fits(outfile,_map_NAIVE_,planckType<double>());
  cout<<"Naive scan map "<<outfile<<endl;

  outfile="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/naive_hits.fits";
  write_Healpix_map_to_fits(outfile,_map_LVF_scan_,planckType<double>());
  cout<<"Naive hits map "<<outfile<<endl;
}

void postprocessing(string file)
{

  hpint64 nside=512;
  load_fits_map_cxx("/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/naive_scan.fits",_map_NAIVE_);
  load_fits_map_cxx("/home/cfeng/computing/spherex/LVFscans/output/wmap_band_imap_r9_9yr_W_v5.fits",_map_GAL_);

  get_naive_map_powspec();
}




