#include "lvf_scans.h"

#include "global.h"

#include "systematics.h"

void copy_vec2vec(double *vec_object, double *vec_object_new)
{
  int k;
  for(k=0;k<3;k++) vec_object_new[k]=vec_object[k];
}

void copy_vec2mat(double *vec_object, double **mat_object)
{
  int k;
  for(k=0;k<3;k++) mat_object[k][0]=vec_object[k];
}
void copy_mat2vec(double **mat_object, double *vec_object)
{
  int k;
  for(k=0;k<3;k++) vec_object[k]=mat_object[k][0];
}

void move_object(double **mat, double *pos, double *object, double *object_new)
{
  int k;
  double **vec3=create_2d_grid<double>(3,1);
  double **vec3_new=create_2d_grid<double>(3,1);

  copy_vec2mat(object,vec3);
  AdotB(mat,vec3,vec3_new,3,3,1);

  for(k=0;k<3;k++) vec3_new[k][0]+=pos[k];

  copy_mat2vec(vec3_new,object_new);

  delete_2d_grid<double>(vec3,3,1);
  delete_2d_grid<double>(vec3_new,3,1);
}

void instrument(double *center_LVF, double **mat_SKY_pointing_0, pair_mat *pm)
{
  double pix=dtr->pix;
  double bottom_vertex[6][3];
  double top_vertex[6][3];
  double bottom_vertex_offset[6][3];
  int num_vertex=6;

  int i,j,k;
  double rb,rt;
  double height,height_box,height_cylinder;

  double r_det;
  int num_det;
  double dphi;
  double phi;


  double p1[3],p2[3],p3[3],p4[3];
  int ii,jj;

  double Fac_larger=1.0;

  rb=5.0*deg2rad*Fac_larger;
  rt=15.0*deg2rad*Fac_larger;
  height=6.0*deg2rad*Fac_larger;
  height_box=1.0*deg2rad*Fac_larger;
  r_det=2.0*deg2rad*Fac_larger;
  num_det=50;
  height_cylinder=10.0*deg2rad*Fac_larger;

  double **inner_cylinder_t=create_2d_grid<double>(num_det,3);
  double **inner_cylinder_b=create_2d_grid<double>(num_det,3);
  double **inner_cylinder_t1=create_2d_grid<double>(num_det,3);
  double **inner_cylinder_b1=create_2d_grid<double>(num_det,3);

  double **mat_y=create_2d_grid<double>(3,3);
  double **vec_cylinder_grid=create_2d_grid<double>(3,1);
  double **vec_cylinder_grid_new=create_2d_grid<double>(3,1);
  double vec_zero[3];//=new create_2d_grid<double>(3,1);


  for(i=0;i<num_vertex;i++)
  {
    phi=60*i;
    phi*=deg2rad;

    bottom_vertex[i][0]=rb*cos(phi);
    bottom_vertex[i][1]=rb*sin(phi);
    bottom_vertex[i][2]=0.0;

    bottom_vertex_offset[i][0]=bottom_vertex[i][0];
    bottom_vertex_offset[i][1]=bottom_vertex[i][1];
    bottom_vertex_offset[i][2]=-height_box;

    top_vertex[i][0]=rt*cos(phi);
    top_vertex[i][1]=rt*sin(phi);
    top_vertex[i][2]=height;
  }

//get the boresight

  dphi=360.0/num_det;



  for(k=0;k<3;k++) vec_zero[k]=0.0;


  double cant_angle=TEST_CANT;

  rot_x(cant_angle*deg2rad,mat_y);


  double **mat_SKY_pointing=create_2d_grid<double>(3,3);
  double **tmp=create_2d_grid<double>(3,3);
  AdotB(pm->mat_phi,pm->mat_theta,tmp,3,3,3);
  AdotB(tmp,mat_SKY_pointing_0,mat_SKY_pointing,3,3,3);

  double boresight_b[3],boresight_t[3];


#ifdef _USE_OPENGL_
  for(i=0;i<num_det;i++)
  {
    phi=i*dphi;

    inner_cylinder_b[i][0]=r_det*cos(phi);
    inner_cylinder_b[i][1]=r_det*sin(phi);
    inner_cylinder_b[i][2]=0;

    inner_cylinder_t[i][0]=inner_cylinder_b[i][0];
    inner_cylinder_t[i][1]=inner_cylinder_b[i][1];
    inner_cylinder_t[i][2]=height_cylinder;


    //double *p=inner_cylinder_b[i];
    move_object(mat_y,vec_zero,inner_cylinder_b[i],inner_cylinder_b[i]);
    move_object(mat_y,vec_zero,inner_cylinder_t[i],inner_cylinder_t[i]);

  }
  
  boresight_b[0]=0;
  boresight_b[1]=0;
  boresight_b[2]=0;

  boresight_t[0]=0;
  boresight_t[1]=0;
  boresight_t[2]=1.5*height_cylinder;
  move_object(mat_y,vec_zero,boresight_b,boresight_b);
  move_object(mat_y,vec_zero,boresight_t,boresight_t);




/*
bottom_vertex
bottom_vertex_offset
top_vertex
inner_cylinder_t
inner_cylinder_b
*/

  for(i=0;i<num_vertex;i++)
  {
    move_object(mat_SKY_pointing,center_LVF,bottom_vertex[i],bottom_vertex[i]);
    move_object(mat_SKY_pointing,center_LVF,bottom_vertex_offset[i],bottom_vertex_offset[i]);
    move_object(mat_SKY_pointing,center_LVF,top_vertex[i],top_vertex[i]);
  }


  for(i=0;i<num_det;i++)
  {
    move_object(mat_SKY_pointing,center_LVF,inner_cylinder_b[i],inner_cylinder_b[i]);
    move_object(mat_SKY_pointing,center_LVF,inner_cylinder_t[i],inner_cylinder_t[i]);
  }
  move_object(mat_SKY_pointing,center_LVF,boresight_b,boresight_b);
  move_object(mat_SKY_pointing,center_LVF,boresight_t,boresight_t);
#endif

  double height_LVF=FAR_VIEW*1.07;
  double LVF_center[3];

  double FOV_theta_0=Shoot_LVF(mat_SKY_pointing, mat_y, center_LVF, height_LVF, LVF_center);

#ifdef _USE_OPENGL_
/*
  cout<<"Instrument->"<<LVF_center[0]<<"\t"<<LVF_center[1]<<"\t"<<LVF_center[2]<<endl;

  get_LVF_center_circles(LVF_center);

  if(FOV_theta_0<=POLE_BOUND) 
  {
    //SHOW_detector2boresight(boresight_t);

    make_movie(Count_NCP_frames,"scan/NCP",windW,windH);

    Count_NCP_frames++;
  }

*/
#endif
  




//start the visualization




#ifdef _USE_OPENGL_
  reset_color();

  glPushMatrix();
  
  for(i=0;i<num_vertex;i++)
  {
    glColor3f(1,1.0/(i+1.0),0);

    if(i==(num_vertex-1)) ii=0;
    else ii=i+1;

    copy_vec2vec(bottom_vertex[i],p1);
    copy_vec2vec(bottom_vertex[ii],p2);
    copy_vec2vec(top_vertex[ii],p3);
    copy_vec2vec(top_vertex[i],p4);


    glBegin(GL_QUADS);
      glVertex3f(p1[0],p1[1],p1[2]);
      glVertex3f(p2[0],p2[1],p2[2]);
      glVertex3f(p3[0],p3[1],p3[2]);
      glVertex3f(p4[0],p4[1],p4[2]);
    glEnd();
  }
  glPopMatrix();


  glPushMatrix();
  glColor3f(0,1,1);
  for(i=0;i<num_vertex;i++)
  {
    if(i==(num_vertex-1)) ii=0;
    else ii=i+1;

    copy_vec2vec(bottom_vertex[i],p1);
    copy_vec2vec(bottom_vertex[ii],p2);
    copy_vec2vec(bottom_vertex_offset[ii],p3);
    copy_vec2vec(bottom_vertex_offset[i],p4);

    glBegin(GL_QUADS);
      glVertex3f(p1[0],p1[1],p1[2]);
      glVertex3f(p2[0],p2[1],p2[2]);
      glVertex3f(p3[0],p3[1],p3[2]);
      glVertex3f(p4[0],p4[1],p4[2]);
    glEnd();
  }
  glPopMatrix();

  glPushMatrix();
  glColor3f(0,1,0.2);
  for(i=0;i<num_det;i++)
  {
    if(i==(num_det-1)) ii=0;
    else ii=i+1;

    copy_vec2vec(inner_cylinder_b[i],p1);
    copy_vec2vec(inner_cylinder_b[ii],p2);
    copy_vec2vec(inner_cylinder_t[ii],p3);
    copy_vec2vec(inner_cylinder_t[i],p4);

    glBegin(GL_QUADS);
      glVertex3f(p1[0],p1[1],p1[2]);
      glVertex3f(p2[0],p2[1],p2[2]);
      glVertex3f(p3[0],p3[1],p3[2]);
      glVertex3f(p4[0],p4[1],p4[2]);
    glEnd();
  }
  glPopMatrix();

  glPushMatrix();
  glColor3f(1,0,0);
  glBegin(GL_LINES);
    glVertex3f(boresight_b[0],boresight_b[1],boresight_b[2]);
    glVertex3f(boresight_t[0],boresight_t[1],boresight_t[2]);
  glEnd();
  glPopMatrix();
#endif



  delete_2d_grid<double>(inner_cylinder_t,num_det,3);
  delete_2d_grid<double>(inner_cylinder_b,num_det,3);
  delete_2d_grid<double>(inner_cylinder_t1,num_det,3);
  delete_2d_grid<double>(inner_cylinder_b1,num_det,3);

  delete_2d_grid<double>(mat_y,3,3);
  delete_2d_grid<double>(vec_cylinder_grid,3,1);
  delete_2d_grid<double>(vec_cylinder_grid_new,3,1);


  delete_2d_grid<double>(tmp,3,3);
  delete_2d_grid<double>(mat_SKY_pointing,3,3);




}

double Shoot_LVF(double **mat_inst, double **mat_boresight, double *pos, double height_LVF, double *boresight_t_OUT)
{
  int i,j,k;
  int ii,jj;

  double strip_width;

  double signal, noise;
 
  double FOV_theta,FOV_phi;

  int numx,numy;
  double xmin,ymin,pix;
  double x,y,z;
  double FOV[3],FOV_sky[3],FOV_sky_unit[3];
  double vec_zero[3];

  double Enlarge_fac=FAR_VIEW;

  for(k=0;k<3;k++) vec_zero[k]=0.0;

  double dis;
  double view_center_ra;
  double view_center_dec;

  hpint64 nside,idpix,id_parent;

  nside=_nside_;

  numx=dtr->numx;
  numy=dtr->numy;
  xmin=dtr->xmin;
  ymin=dtr->ymin;
  pix=dtr->pix;
  strip_width=dtr->strip_width;

  double FOV_theta_0,FOV_phi_0;
  x=xmin+pix*numx/2;
  y=ymin+pix*numy/2;
  z=height_LVF;
  FOV[0]=x*Enlarge_fac;
  FOV[1]=y*Enlarge_fac;
  FOV[2]=z;
  move_object(mat_boresight,vec_zero,FOV,FOV);
  move_object(mat_inst,pos,FOV,FOV);

  dis=get_vec_len(FOV);
  if(dis!=0) for(i=0;i<3;i++) FOV_sky_unit[i]=FOV[i]/dis;
  else for(i=0;i<3;i++) FOV_sky_unit[i]=1e-8;
  spherical_2_cartesian(FOV_sky_unit[0],FOV_sky_unit[1],FOV_sky_unit[2],FOV_theta_0,FOV_phi_0);


if(_ID_STEP_%8==0)
//if(FOV_theta_0<=(1.5*POLE_BOUND))
{
/*
  for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
  {
    x=xmin+pix*ii;
    y=ymin+pix*jj;
    z=height_LVF;

    FOV[0]=x*Enlarge_fac;
    FOV[1]=y*Enlarge_fac;
    FOV[2]=z;

    move_object(mat_boresight,vec_zero,FOV,FOV);
    move_object(mat_inst,pos,FOV,FOV);



    double dis=get_vec_len(FOV);
    if(dis!=0) for(i=0;i<3;i++) FOV_sky_unit[i]=FOV[i]/dis;
    else for(i=0;i<3;i++) FOV_sky_unit[i]=1e-8;


    int idband=y/strip_width;
    RECORD_NCP_flat(idband,FOV_sky_unit);

    spherical_2_cartesian(FOV_sky_unit[0],FOV_sky_unit[1],FOV_sky_unit[2],FOV_theta,FOV_phi);

    ang2pix_ring64(nside,FOV_theta*deg2rad,FOV_phi*deg2rad,&idpix);

    signal=_map_GAL_[idpix];
    double eta=_rng_.rand_gauss();
    double noise_level=2.0;
    noise=eta*noise_level;

    dtr->fov_coord[ii][jj].ra=x;
    dtr->fov_coord[ii][jj].dec=y;
    dtr->fov_image[ii][jj]=signal;

    dtr->fov_coord3d_x[ii][jj]=FOV[0];
    dtr->fov_coord3d_y[ii][jj]=FOV[1];
    dtr->fov_coord3d_z[ii][jj]=FOV[2];


    _map_NAIVE_[idpix]+=(signal+noise);
    _map_LVF_scan_[idpix]++;

    if( (ii==(numx/2)) && (jj==(numy/2)) ) 
    {
      view_center_ra=FOV_phi;
      view_center_dec=FOV_theta;

      if(view_center_dec<=POLE_BOUND) Count_NCP_pointings++;
#ifdef _USE_OPENGL_
      cout<<"View center -> ("<<view_center_ra<<","<<view_center_dec<<")"<<" Count_NCP_pointings->"<<Count_NCP_pointings<<endl;
#endif
    }
*/
 

    Count_NCP_frames++;
 

cout<<"DEEP day: "<<Count_NCP_frames<<"\t"<<_ID_STEP_<<endl;
  //for(k=0;k<96;k++)
  for(k=0;k<144;k++)
  {
//if((k==11)||(k==35)) 
//if((k>=40)&&(k<=56)) continue;
    //if( ((k>=24)&&(k<=71)) || ((k>=102)&& (k<=137)) ) continue;

//    if( ((k>=12)&&(k<=35))||((k>=60)&&(k<=83)) )
//if(k==83)

   {
    shoot_deep(k,mat_inst,DEEP_pointings[k][0],DEEP_pointings[k][1]-PI/2.0);
    //cout<<"DEEP:"<<k<<endl;
  
//cout<<ii<<"\t"<<jj<<"\t"<<FOV[0]<<"\t"<<
  






#ifdef _USE_OPENGL_
  double boresight_b[3],boresight_t[3];
  boresight_b[0]=0.0;
  boresight_b[1]=0.0;
  boresight_b[2]=0.0;

  boresight_t[0]=0.0;
  boresight_t[1]=0.0;
  boresight_t[2]=height_LVF;

  move_object(mat_boresight,vec_zero,boresight_b,boresight_b);
  move_object(mat_inst,pos,boresight_b,boresight_b);

  move_object(mat_boresight,vec_zero,boresight_t,boresight_t);
  move_object(mat_inst,pos,boresight_t,boresight_t);

  copy_vec2vec(boresight_t,boresight_t_OUT);


  double rr,gg,bb;
  double v;
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

  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex3f(boresight_b[0],boresight_b[1],boresight_b[2]);
  glVertex3f(boresight_t[0],boresight_t[1],boresight_t[2]);
  glEnd();
  glPopMatrix();
#endif
}
}

if(Count_NCP_frames%30==0)
{
  Analyze_NCP_flat(Count_NCP_frames);
  Analyze_NCP_image(Count_NCP_frames);


  char seg[1024];
  memset(seg,0,1024);
  sprintf(seg,"%03d",Count_NCP_frames);
  string file,null="";
  string path="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";
  file=null+"LeastHits_"+seg;
  //get_angular_counts(path,file);
}
  //OUTPUT_NCP_flat(Count_NCP_frames);
//make_movie(Count_NCP_frames,"scan/NCP",windW,windH);


}

  return FOV_theta_0;

}


int get_LVF_band(const double x, const double y)
{
  int numx,numy;
  int id;
  double xmin,ymin,pix,strip_width;
  numx=dtr->numx;
  numy=dtr->numy;
  xmin=dtr->xmin;
  ymin=dtr->ymin;
  pix=dtr->pix;
  strip_width=dtr->strip_width;

  if(y<0)
  {
    id=(x-xmin)/strip_width;
  }
  else
  {
    id=(x-xmin)/strip_width;
    id=47-id;
  }

//  id=(y-ymin)/strip_width;

  return id;
}


//theta, phi [rad]
void shoot_deep(const int id_pts,double **Amat, const double theta, const double phi)
{
  //double **mat_y=create_2d_grid<double>(3,3);
  //double **mat_z=create_2d_grid<double>(3,3);
  //double **Amat=create_2d_grid<double>(3,3);

  double signal, noise; 
  double FOV_theta,FOV_phi;

  int numx,numy;
  double xmin,ymin,pix;
  double x,y,z;
  double FOV[3],FOV_sky[3],FOV_sky_unit[3];
  double vec_zero[3];
  int ii,jj,i,j,k;
  double strip_width;

  hpint64 nside,idpix,id_parent;

  nside=_nside_;

  numx=dtr->numx;
  numy=dtr->numy;
  xmin=dtr->xmin;
  ymin=dtr->ymin;
  pix=dtr->pix;
  strip_width=dtr->strip_width;

  double Enlarge_fac=FAR_VIEW;
  double height_LVF=FAR_VIEW*1.07;

/*
  rot_y(theta,mat_y);
  rot_z(phi,mat_z);
  AdotB(mat_z,mat_y,Amat,3,3,3);
*/
  double pos[3];//=new double [3];
  double zero[3];

  double r=Enlarge_fac;

  zero[0]=0;
  zero[1]=0;
  zero[2]=0;

  pos[0]=r*sin(theta)*cos(phi);
  pos[1]=r*sin(theta)*sin(phi);
  pos[2]=r*cos(theta);

//#ifdef _rotate_lvf_
//  ROTATE_LVF_90(pos);

//#endif
  double pos_running[3];
  move_object(Amat,zero,pos,pos_running);


//ofstream out("COLORS");

  for(ii=0;ii<numx;ii++) for(jj=0;jj<numy;jj++)
  {
    x=xmin+pix*ii;
    y=ymin+pix*jj;
    z=0;
#ifdef _rotate_lvf_
double aa[3];
aa[0]=x;
aa[1]=y;
aa[2]=z;
ROTATE_LVF_90(aa);
x=aa[0];
y=aa[1];
z=aa[2];
#endif
    FOV[0]=x*Enlarge_fac;
    FOV[1]=y*Enlarge_fac;
    FOV[2]=z;


    move_object(Amat,pos_running,FOV,FOV);

    double dis=get_vec_len(FOV);
    if(dis!=0) for(i=0;i<3;i++) FOV_sky_unit[i]=FOV[i]/dis;
    else for(i=0;i<3;i++) FOV_sky_unit[i]=1e-8;

    int idband;
    //int idband=(y-ymin)/strip_width;    

    idband=get_LVF_band(x,y);

//out<<x<<"\t"<<y<<"\t"<<idband<<endl;
/*
#ifdef _rotate_lvf_
    idband=(x-ymin)/strip_width;
#else
    idband=(y-ymin)/strip_width;
#endif
*/

    RECORD_NCP_flat(idband,FOV_sky_unit);

    spherical_2_cartesian(FOV_sky_unit[0],FOV_sky_unit[1],FOV_sky_unit[2],FOV_theta,FOV_phi);

    ang2pix_ring64(nside,FOV_theta*deg2rad,FOV_phi*deg2rad,&idpix);

    signal=_map_GAL_[idpix];
    double eta=_rng_.rand_gauss();
    double omega=(6.2/3600.0*deg2rad); omega=omega*omega;
    double noise_level=1.35;//sqrt(2e-8/omega);
    noise=eta*noise_level;


    
    //double ZL=get_ZL(id_pts,x,y);
    double MEVconversion=414.0;
    double dark_current=MEVconversion*mat_dark_current_LOW[ii][jj];//get_stamp_pixel(x,y,mat_dark_current_LOW);

    double rnd_thermo=_rng_.rand_gauss();
    double dark_thermo=MEVconversion*mat_thermo_LOW[ii][jj];//get_stamp_pixel(x,y,mat_thermo_LOW);
    if(dark_thermo<0) dark_thermo=0;

    dark_thermo*=rnd_thermo;

    double zodi;
    int id_step=id_pts%4;
    switch(id_step)
    {
      case 0:
        zodi=mat_zodi_A_LOW[ii][jj];
        break;
      case 1:
        zodi=0.0;
        break;
      case 2:
        zodi=0.0;
        break;
      case 3:
        zodi=mat_zodi_B_LOW[ii][jj];
        break;
    }


    double all_components=signal;//414*(dark_thermo);

//broad band 1.6 with R=4
if((idband>=24)&&(idband<35))
{
    capture_NCP_image(0,signal,FOV_sky_unit);
    capture_NCP_image(1,noise,FOV_sky_unit);
    capture_NCP_image(2,zodi,FOV_sky_unit);
    capture_NCP_image(3,dark_current,FOV_sky_unit);
    capture_NCP_image(4,dark_thermo,FOV_sky_unit);
}


    dtr->fov_coord[ii][jj].ra=x;
    dtr->fov_coord[ii][jj].dec=y;
    dtr->fov_image[ii][jj]=all_components;


    dtr->fov_coord3d_x[ii][jj]=FOV[0];
    dtr->fov_coord3d_y[ii][jj]=FOV[1];
    dtr->fov_coord3d_z[ii][jj]=FOV[2];
  }

//out.close();

//throw;
  //delete_2d_grid<double>(mat_y,3,3);
  //delete_2d_grid<double>(mat_z,3,3);
  //delete_2d_grid<double>(Amat,3,3);
  //delete pos;
}

void Shoot_LVF_fast(double **mat_inst, double **mat_boresight, double *pos, double height_LVF, double *boresight_t_OUT)
{
  int i,j,k;
  int ii,jj;

  double signal, noise, eta;
  double noise_level=2.0;
 
  double FOV_theta,FOV_phi;

  int numx,numy;
  double xmin,ymin,pix;
  double x,y,z;
  double FOV[3],FOV_sky[3],FOV_sky_unit[3];
  double vec_zero[3];

  double Enlarge_fac=FAR_VIEW;

  for(k=0;k<3;k++) vec_zero[k]=0.0;

  double dis;
  double view_center_ra;
  double view_center_dec;
  double theta, phi;

  hpint64 nside,idpix,id_parent;

  nside=_nside_;

  numx=dtr->numx;
  numy=dtr->numy;
  xmin=dtr->xmin;
  ymin=dtr->ymin;
  pix=dtr->pix;

  x=xmin+pix*numx/2;
  y=ymin+pix*numy/2;
  z=height_LVF;



  FOV[0]=x*Enlarge_fac;
  FOV[1]=y*Enlarge_fac;
  FOV[2]=z;


  move_object(mat_boresight,vec_zero,FOV,FOV);
  move_object(mat_inst,pos,FOV,FOV);

  dis=get_vec_len(FOV);
  if(dis!=0) for(i=0;i<3;i++) FOV_sky_unit[i]=FOV[i]/dis;
  else for(i=0;i<3;i++) FOV_sky_unit[i]=1e-8;
  spherical_2_cartesian(FOV_sky_unit[0],FOV_sky_unit[1],FOV_sky_unit[2],FOV_theta,FOV_phi);
    //ang2pix_ring64(nside,FOV_theta*deg2rad,FOV_phi*deg2rad,&idpix);


  int flag=0;
  double view_radius;
  if(FOV_theta<=POLE_BOUND) 
  {
    Count_NCP_pointings++;


  view_radius=2.5*deg2rad;


  pointing pt;
  pt.theta=PI/2.0-FOV_theta;
  pt.phi=FOV_phi;
  vector<int> listpix;
  _map_GAL_.query_disc_inclusive(pt,view_radius,listpix);




  //cout<<"Shoot LVF fast->"<<listpix.size()<<"\t"<<FOV_theta<<"\t"<<FOV_phi<<endl;
  for(i=0;i<listpix.size();i++)
  {
    idpix=listpix[i];

    signal=_map_GAL_[idpix];
    eta=_rng_.rand_gauss();
    noise=eta*noise_level;

    _map_NAIVE_[idpix]+=(signal+noise);
    _map_LVF_scan_[idpix]++;
    
    //pix2ang_ring64(nside, idpix, &theta, &phi);


    //dtr->fov_coord[ii][jj].ra=phi;
    //dtr->fov_coord[ii][jj].dec=theta;
    //dtr->fov_image[ii][jj]=signal;
/*
    dis=FAR_VIEW;
    FOV[0]=dis*sin(theta)*cos(phi);
    FOV[1]=dis*sin(theta)*cos(phi);
    FOV[2]=dis*cos(theta);

    dtr->fov_coord3d_x[ii][jj]=FOV[0];
    dtr->fov_coord3d_y[ii][jj]=FOV[1];
    dtr->fov_coord3d_z[ii][jj]=FOV[2];
*/

  }

}
  //if(FOV_theta<=POLE_BOUND) Count_NCP_pointings++;




#ifdef _USE_OPENGL_
  double boresight_b[3],boresight_t[3];
  boresight_b[0]=0.0;
  boresight_b[1]=0.0;
  boresight_b[2]=0.0;

  boresight_t[0]=0.0;
  boresight_t[1]=0.0;
  boresight_t[2]=height_LVF;

  move_object(mat_boresight,vec_zero,boresight_b,boresight_b);
  move_object(mat_inst,pos,boresight_b,boresight_b);

  move_object(mat_boresight,vec_zero,boresight_t,boresight_t);
  move_object(mat_inst,pos,boresight_t,boresight_t);

  copy_vec2vec(boresight_t,boresight_t_OUT);


  double rr,gg,bb;
  double v;
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

  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex3f(boresight_b[0],boresight_b[1],boresight_b[2]);
  glVertex3f(boresight_t[0],boresight_t[1],boresight_t[2]);
  glEnd();
  glPopMatrix();
#endif

}

void SHOW_detector2boresight(double *pos)
{
  double corner[4][3];
  double center_LVF[3];
  int i,j,k;
  int nx=dtr->numx;
  int ny=dtr->numy;


  copy_vec2vec(pos,center_LVF);

  corner[0][0]=dtr->fov_coord3d_x[0][0];
  corner[0][1]=dtr->fov_coord3d_y[0][0];
  corner[0][2]=dtr->fov_coord3d_z[0][0];

  corner[1][0]=dtr->fov_coord3d_x[nx-1][0];
  corner[1][1]=dtr->fov_coord3d_y[nx-1][0];
  corner[1][2]=dtr->fov_coord3d_z[nx-1][0];

  corner[2][0]=dtr->fov_coord3d_x[nx-1][ny-1];
  corner[2][1]=dtr->fov_coord3d_y[nx-1][ny-1];
  corner[2][2]=dtr->fov_coord3d_z[nx-1][ny-1];

  corner[3][0]=dtr->fov_coord3d_x[0][ny-1];
  corner[3][1]=dtr->fov_coord3d_y[0][ny-1];
  corner[3][2]=dtr->fov_coord3d_z[0][ny-1];


  double p1[3],p2[3],p3[3],p4[3];
  double alpha=0.4;
  reset_color();
  for(k=0;k<3;k++) 
  {
    p1[k]=corner[0][k];
    p2[k]=corner[1][k];
    p3[k]=corner[2][k];
    p4[k]=corner[3][k];

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


struct_lvf_track * new_LVF_track(const int num_step, const int num_grid)
{
  int i;
  struct_lvf_track *p;
  p=new struct_lvf_track;
  p->num_step=num_step;
  p->num_grid=num_grid;
  p->STEP_RING=new COORD3D * [num_step];
  for(i=0;i<num_step;i++) p->STEP_RING[i]=new COORD3D [num_grid];

  p->current_num_rings=0;
  cout<<"Structure LVF_track is allocated with dims ("<<num_step<<","<<num_grid<<")"<<endl;
  return p;
}

void delete_LVF_track(struct_lvf_track *p, const int num_step, const int num_grid)
{
  int i;
  for(i=0;i<num_step;i++) delete p->STEP_RING[i];
  delete p;
}

void insert_LVF_track(COORD3D *ringh, COORD3D *ringv)
{
  int i,j,k;
  int num_grid=LVF_track->num_grid;
  int id_pointing=_ID_STEP_/NUM_STEPPING;
  int id_step=_ID_STEP_%NUM_STEPPING;

  int max_count=2*NUM_STEPPING;

  //if(LVF_track->current_num_rings==max_count) LVF_track->current_num_rings=0;
  if(id_step==0) LVF_track->current_num_rings=0;
//cout<<num_grid<<endl;
  for(i=0;i<num_grid;i++) for(k=0;k<3;k++) LVF_track->STEP_RING[LVF_track->current_num_rings][i].r[k]=ringh[i].r[k];
  LVF_track->current_num_rings++;

  for(i=0;i<num_grid;i++) for(k=0;k<3;k++) LVF_track->STEP_RING[LVF_track->current_num_rings][i].r[k]=ringv[i].r[k];
  LVF_track->current_num_rings++;

#ifdef _USE_OPENGL_
  cout<<"insert_LVF_track->NUM rings: "<<LVF_track->current_num_rings<<endl;
#endif


}

void get_LVF_center_circles(double *pos)
{
  int i,j,k;
  double dis;
  double unit_pos[3];
  double FOV_theta,FOV_phi;
  dis=get_vec_len(pos);
  if(dis!=0) for(i=0;i<3;i++) unit_pos[i]=pos[i]/dis;
  else for(i=0;i<3;i++) unit_pos[i]=1e-8;
  spherical_2_cartesian(unit_pos[0],unit_pos[1],unit_pos[2],FOV_theta,FOV_phi);
#ifdef _USE_OPENGL_
  cout<<"get_LVF_center_circles->"<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<endl;
  cout<<"get_LVF_center_circles->"<<FOV_theta<<"\t"<<FOV_phi<<endl;
#endif

  double radius=dis;
  int num_grid=LVF_track->num_grid;
  COORD3D *ring_v=new COORD3D [num_grid];
  COORD3D *ring_h=new COORD3D [num_grid];

  make_circle_h(radius, FOV_theta, num_grid, ring_h);
  make_circle_v(radius, FOV_phi, num_grid, ring_v);

  insert_LVF_track(ring_h,ring_v);



#ifdef _USE_OPENGL_
  glPushMatrix();
  for(k=0;k<LVF_track->current_num_rings;k++)
  {
    glBegin(GL_LINES);
    for(i=0;i<num_grid;i++)
    {
      glColor3f(1.0,1.0,1.0);
      glVertex3f(LVF_track->STEP_RING[k][i].r[0],LVF_track->STEP_RING[k][i].r[1],LVF_track->STEP_RING[k][i].r[2]);
    }
    glEnd();
  }
  glPopMatrix();

#endif


  delete ring_v;
  delete ring_h;
}


void ROTATE_LVF_90(double *pos)
{
  double pos_copy[3];
  int i;
  for(i=0;i<3;i++) pos_copy[i]=pos[i];

  double **mat_z=create_2d_grid<double>(3,3);
  double **mat_pos1=create_2d_grid<double>(3,1);
  double **mat_pos2=create_2d_grid<double>(3,1);

  for(i=0;i<3;i++) mat_pos1[i][0]=pos[i];

  rot_z(90.0*deg2rad,mat_z);
  AdotB(mat_z,mat_pos1,mat_pos2,3,3,1);

  for(i=0;i<3;i++) pos[i]=mat_pos2[i][0];

}
