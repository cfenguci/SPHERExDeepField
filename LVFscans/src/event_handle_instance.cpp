#include "lvf_scans.h"

#include "global.h"


void get_position(const int out_fits_file)
{
  int max_days=365;
  if(_cnt_days_>10*max_days) return;

  //cout<<"(get_position)Camera position: "<<CAMERA_X<<"\t"<<CAMERA_Y<<"\t"<<CAMERA_Z<<endl;
  
  int i,j,k;
  int ii,jj;
  double t_sky,t_LVF,tmin,tmax,dt;
  int num;
  double x1[3],x2[3],xdiff[3];
  double xn[3],xlong[3];

  double x1_LVF[3],xn_LVF[3];
  double rot[3][3];

  tmin=0.0;
  //tmax=24.0;

  //dt=95.0/3600.0;//185/3600.0;
  dt=12.0/60.0;
  //dt=185.0/3600.0;


  num=_num_orbit_*24/dt;//(tmax-tmin)/dt;


  double height=5e5;
  double v=sqrt(_G_*M_earth/(R_earth+height));
  double T=96*60;//2.0*PI*(R_earth+height)/v;
  double omega_theta=360.0/(T/3600.0);
  double omega_phi;
//cout<<T<<endl;
#ifdef _USE_OPENGL_
cout<<"Orbit time->"<<T<<endl;
  omega_phi=1.0/(96.0/60)*10;
#else
  //omega_phi=1.0/24.0;
  omega_phi=1.0/(96.0/60);
#endif;


  double w=WIDTH_DETECTOR*deg2rad;
  double h=HEIGHT_DETECTOR*deg2rad;
  
  double strip_width=dtr->strip_width;


  double xmin=-h/2.0,xmax=h/2.0;
  double ymin,ymax;
  get_strip(CURRENT_STRIP_ID,NUM_STRIP,ymin,ymax);
  
  strip_width/=deg2rad;
  


  double x,y,z;
  double pix=PIX_REDUC*6.2/60.0/60.0*deg2rad;

  int numx,numy;
  numx=(xmax-xmin)/pix;
  numy=(ymax-ymin)/pix;


  double **mat_x=create_2d_grid<double>(3,3);
  double **mat_y=create_2d_grid<double>(3,3);
  double **mat_z=create_2d_grid<double>(3,3);
  double **Amat=create_2d_grid<double>(3,3);
  double **Amat_LVF=create_2d_grid<double>(3,3);
  double **Bmat=create_2d_grid<double>(3,1);
  double **Cmat=create_2d_grid<double>(3,1);

  
  vec3 vec;
  
  

#ifdef _USE_OPENGL_

  reset_color();
  //glRotatef(-90,1,0,0);
  //glRotatef(-90,0,0,1);
#endif


  
  double theta_out, phi_out;
  double theta_out_LVF, phi_out_LVF;




  int id_slews;
  i=_ID_STEP_;
//cout<<"STEP:"<<i<<endl;
  
  int id=i/NUM_STEPPING;
  id_slews=i%NUM_STEPPING;

  double theta_out_slews, phi_out_slews;


  int RecordSky;

  t_sky=tmin+NUM_STEPPING*dt*id;
  t_LVF=tmin+dt*i;


  double theta_start=INIT_THETA;
  double phi_start=INIT_PHI;
  generate_pointing(1.07, theta_start, phi_start, omega_theta,omega_phi,t_sky,x1,xn,theta_out,phi_out);
  generate_pointing(1.07, theta_start, phi_start, omega_theta,omega_phi,t_LVF,x1_LVF,xn_LVF,theta_out_LVF,phi_out_LVF);

  dtr->omega_theta=omega_theta;
  dtr->omega_phi=omega_phi;
  dtr->delta_t=dt;


//SCP
  //make_circle_h(1.02, 180-POLE_BOUND, 360, 1,0,0);
  
//NEP
  //make_circle_h(1.02, POLE_BOUND, 360, 1,0,0);
  //make_circle_h(3*1.02, POLE_BOUND, 360, 1,0,0);
  


  make_circle(1.07, phi_out, 200, 1,1,0);



  rot_y(theta_out_LVF*deg2rad,mat_y);
  rot_z(phi_out_LVF*deg2rad,mat_z);
  AdotB(mat_z,mat_y,Amat_LVF,3,3,3);



  generate_pointing_slews(1.07,theta_out, phi_out,0.0,omega_theta,omega_phi,id_slews,strip_width, x1, theta_out_slews, phi_out_slews);

  rot_y(theta_out_slews*deg2rad,mat_y);
  rot_z(phi_out_slews*deg2rad,mat_z);
  AdotB(mat_z,mat_y,Amat,3,3,3);


#ifdef _USE_OPENGL_
  RecordSky=1;
#else
  if(theta_out_slews<POLE_BOUND) RecordSky=1;//||(theta_out_slews>(180-POLE_BOUND))) RecordSky=1;
  else RecordSky=0;

  //if(RecordSky==0) cout<<"SKIP: "<<_ID_STEP_<<endl;
#endif

  RecordSky=1;

  FOV_image(RecordSky,dtr,Amat_LVF,x1_LVF);


#ifndef  _USE_OPENGL_
/*
 int month=30*num;

  if(i%month==0) 
  {

    cout<<"------------------- Day "<<_cnt_days_<<"-------------------"<<endl;


//    int max_days=365;
    int dd;
//    if(_cnt_days_>max_days) return;

    Healpix_Map<double> daily_map(_nside_,RING,SET_NSIDE);
	
    for(dd=0;dd<daily_map.Npix();dd++)
    {
      double nhits=_map_LVF_scan_[dd];
      if(nhits!=0) daily_map[dd]=_map_NAIVE_[dd]/nhits;
      else daily_map[dd]=0.0;
    }

    get_grid(_cnt_days_,daily_map, "grid_image");
    get_grid(_cnt_days_,_map_LVF_scan_, "grid_hits");



    char seg[1024];
    memset(seg,0,1024);
    sprintf(seg,"%d",_cnt_days_);
    string null="";


    get_angular_counts((null+"grid_hits_"+seg+".txt").c_str());
cout<<"Count_NCP_pointings->"<<Count_NCP_pointings<<endl;

    _cnt_days_+=30;
  }

*/
#endif


  if(out_fits_file)
  {
    string outfile="/home/cmb/computing/spherex/VLFscans/output/VLF_scan.fits";
    write_Healpix_map_to_fits(outfile,_map_LVF_scan_,planckType<double>());
  }

  delete_2d_grid<double>(mat_x,3,3);
  delete_2d_grid<double>(mat_y,3,3);
  delete_2d_grid<double>(mat_z,3,3);
  delete_2d_grid<double>(Amat,3,1);
  delete_2d_grid<double>(Amat_LVF,3,1);
  delete_2d_grid<double>(Bmat,3,1);
  delete_2d_grid<double>(Cmat,3,1);

}

