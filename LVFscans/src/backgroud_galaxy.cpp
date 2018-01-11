#include "lvf_scans.h"

#include "global.h"

void test_color()
{
  int i;
  double r,g,b;
  for(i=0;i<256*256*256;i++)
  {
    color2rgb(i,r,g,b);
    cout<<r<<"\t"<<g<<"\t"<<b<<endl;
  }
}

void color2rgb(const int c, double &rr, double &gg, double &bb)
{
  int r=(c>>16)&255;
  int g=(c>>8)&255;
  int b=c&255;

  double fac=255.0;

  rr=r/fac;
  gg=g/fac;
  bb=b/fac;
  //cout<<r/fac<<"\t"<<g/fac<<"\t"<<b/fac<<endl; 
}

void get_map_pix_stat()
{
  int i;
  long npix=_map_GAL_.Npix();

  double max=-1e100;
  double min=1e100;
  double v;
  for(i=0;i<npix;i++)
  {
    v=_map_GAL_[i];
    if(v>=max) max=v;
    if(v<=min) min=v;
  }
  cout<<"Map value range: ["<<min<<"->"<<max<<"]"<<endl;

  int nbin=1000;
  double nstep=(max-min)/nbin;

  int *nbin_count=new int [nbin+1];
  for(i=0;i<nbin+1;i++) nbin_count[i]=0;

  for(i=0;i<npix;i++)
  {
    v=_map_GAL_[i];
    int id=(v-min)/nstep;
    nbin_count[id]++;
  }

  ofstream out;
  out.open("/home/cfeng/computing/spherex/LVFscans/output/hist_map_pix");
  for(i=0;i<nbin;i++)
  {
    double v1,v2;
    v1=min+nstep*i;
    v2=v1+nstep;
    out<<(v1+v2)/2.0<<"\t"<<nbin_count[i]<<endl;
  }
  out.close();

  delete nbin_count;
}



hpint64 get_index_mapping_parent2kid(hpint64 id_parent)
{
  int nside=_nside_;
  double theta,phi;
  double rot=_ID_STEP_*185.0/3600.0*360.0/24.0;
  hpint64 id_rot,id_equ;
  double l,b,ra,dec;
  
  vec3 bk_xyz,bk_xyz_new;
  
  double **mat_z=create_2d_grid<double>(3,3);
  double **mat_pix=create_2d_grid<double>(3,1);  
  double **mat_pix_new=create_2d_grid<double>(3,1);  
  
  //pix2ang_ring64(nside, id_kid, &theta, &phi);  
  bk_xyz=_map_GAL_.pix2vec(id_parent);
  
  mat_pix[0][0]=bk_xyz.x;
  mat_pix[1][0]=bk_xyz.y;
  mat_pix[2][0]=bk_xyz.z;
	
  rot_z(rot*deg2rad,mat_z);
  AdotB(mat_z,mat_pix,mat_pix_new,3,3,1);
	
  bk_xyz_new.x=mat_pix_new[0][0];
  bk_xyz_new.y=mat_pix_new[1][0];
  bk_xyz_new.z=mat_pix_new[2][0];
	
	
  id_rot=_map_GAL_.vec2pix(bk_xyz_new);
	
  pix2ang_ring64(nside, id_rot, &theta, &phi);
//cout<<"--->"<<id_parent<<"\t"<<id_rot<<"\t"<<theta/deg2rad<<"\t"<<phi/deg2rad<<"\t";
  l=phi;
  b=PI/2.0-theta;
  gal2equ(l,b,ra,dec);
//cout<<dec/deg2rad<<"\t"<<ra/deg2rad<<"\t";
	
  theta=PI/2.0-dec;
  phi=ra;
	
  ang2pix_ring64(nside,theta,phi,&id_equ);
//cout<<id_equ<<endl;
  
  delete_2d_grid<double>(mat_pix,3,1);
  delete_2d_grid<double>(mat_pix_new,3,1);
  delete_2d_grid<double>(mat_z,3,3);
  
  return id_equ;
}

hpint64 get_index_mapping_kid2parent(hpint64 id_kid)
{
  int nside=_nside_;
  double theta,phi;

  double rot=_ID_STEP_*185.0/3600.0*360.0/24.0;
  
  hpint64 id_derot,id_gal;
  double l,b,ra,dec;
  
  vec3 bk_xyz,bk_xyz_new;
  
  double **mat_z=create_2d_grid<double>(3,3);
  double **mat_pix=create_2d_grid<double>(3,1);  
  double **mat_pix_new=create_2d_grid<double>(3,1);  
  
  
  pix2ang_ring64(nside, id_kid, &dec, &ra);
//cout<<"<---"<<dec/deg2rad<<"\t"<<ra/deg2rad<<"\t";
  dec=PI/2.0-dec;
  equ2gal(ra,dec,l,b);
	
  theta=PI/2.0-b;
  phi=l;
	
  ang2pix_ring64(nside,theta,phi,&id_gal);
//cout<<id_gal<<"\t"<<theta/deg2rad<<"\t"<<phi/deg2rad<<"\t";
  
  //pix2ang_ring64(nside, id_kid, &theta, &phi);  
  bk_xyz=_map_GAL_.pix2vec(id_gal);
  
  mat_pix[0][0]=bk_xyz.x;
  mat_pix[1][0]=bk_xyz.y;
  mat_pix[2][0]=bk_xyz.z;
	
  rot_z(-rot*deg2rad,mat_z);
  AdotB(mat_z,mat_pix,mat_pix_new,3,3,1);
	
  bk_xyz_new.x=mat_pix_new[0][0];
  bk_xyz_new.y=mat_pix_new[1][0];
  bk_xyz_new.z=mat_pix_new[2][0];
	
	
  id_derot=_map_GAL_.vec2pix(bk_xyz_new);
//cout<<id_derot<<endl;
  
  delete_2d_grid<double>(mat_pix,3,1);
  delete_2d_grid<double>(mat_pix_new,3,1);
  delete_2d_grid<double>(mat_z,3,3);

  
  return id_derot;
}

void rotate_background()
{
  int i;
  hpint64 idpix;
  int nside=_nside_;
  double rot=_ID_STEP_*185.0/3600.0*360.0/240.0;
  
  double theta, phi;
  
  double l,b,ra,dec;
  
  hpint64 id_rot,id_equ;
  
  double **mat_pix=create_2d_grid<double>(3,1);
  double **mat_pix_new=create_2d_grid<double>(3,1);
  double **mat_z=create_2d_grid<double>(3,3);
  double **mat_x=create_2d_grid<double>(3,3);
  double **mat_y=create_2d_grid<double>(3,3);
  double **mat_tmp=create_2d_grid<double>(3,3);
  double **mat_tmp1=create_2d_grid<double>(3,3);
  Healpix_Map<double> _map_GAL_new_(_nside_,RING,SET_NSIDE);  
  vec3 bk_xyz,bk_xyz_new;
  
  
  rot_z(0.0,mat_z);
  //AdotB(mat_z,mat_x,mat_tmp1,3,3,3);

  
  for(i=0;i<_map_GAL_.Npix();i++)
  {

    //pix2ang_ring64(nside, i, &theta, &phi);
    //theta/=deg2rad;
    //phi/=deg2rad;
	

	
	bk_xyz=_map_GAL_.pix2vec(i);
	mat_pix[0][0]=bk_xyz.x;
	mat_pix[1][0]=bk_xyz.y;
	mat_pix[2][0]=bk_xyz.z;
	
	
	AdotB(mat_z,mat_pix,mat_pix_new,3,3,1);
	
	bk_xyz_new.x=mat_pix_new[0][0];
	bk_xyz_new.y=mat_pix_new[1][0];
	bk_xyz_new.z=mat_pix_new[2][0];
	
	
	id_rot=_map_GAL_.vec2pix(bk_xyz_new);
	
	pix2ang_ring64(nside, id_rot, &theta, &phi);
	l=phi;
	b=PI/2.0-theta;
	
	gal2equ(l,b,ra,dec);
	
	theta=PI/2.0-dec;
	phi=ra;
	
	ang2pix_ring64(nside,theta,phi,&id_equ);
	
	_index_mapping_[id_equ]=i;
	
	//id_equ=id_rot;
	
	
	
	_map_GAL_new_[id_equ]=_map_GAL_[i];	
  }
  
  for(i=0;i<_map_GAL_.Npix();i++) _map_GAL_[i]=_map_GAL_new_[i];
  cout<<"Map updated! "<<rot<<endl;


  //_map_GAL_256_.Import_degrade(_map_GAL_);
  //cout<<"Background galaxy image with nside=256 is created!"<<endl;
  
  delete_2d_grid<double>(mat_pix,3,1);
  delete_2d_grid<double>(mat_pix_new,3,1);
  delete_2d_grid<double>(mat_z,3,3);
  delete_2d_grid<double>(mat_x,3,3);
  delete_2d_grid<double>(mat_y,3,3);
  delete_2d_grid<double>(mat_tmp,3,3);
  delete_2d_grid<double>(mat_tmp1,3,3);
}

/*
double *RING_position2vec(const double l, const double b);
{
  double ll,bb;

  ll=l*deg2rad;
  bb=b*deg2rad;
  bb=PI/2.0-bb;

  double ra,dec;
  double *x=new double [3];

  gal2equ(ll,bb,ra,dec);

  dec=PI/2.0-dec;

  x[0]=radius*sin(dec)*cos(ra);
  x[1]=radius*sin(dec)*sin(ra);
  x[2]=radius*cos(dec);

  return x;
}
*/

void make_circle_galactic(const double radius, const double theta, const int num_grid, COORD3D *ring)
{
  int i;
  double ra,dec,l,b;
  double dphi=360.0/num_grid;
  double x[3];
  //glPushMatrix();
  //glBegin(GL_LINES);
  for(i=0;i<num_grid;i++)
  {

    l=i*dphi;
    b=theta;




    l*=deg2rad;
    b*=deg2rad;
    b=PI/2.0-b;

    gal2equ(l,b,ra,dec);

    dec=PI/2.0-dec;

    x[0]=radius*sin(dec)*cos(ra);
    x[1]=radius*sin(dec)*sin(ra);
    x[2]=radius*cos(dec);


    for(int k=0;k<3;k++) ring[i].r[k]=x[k];

    //glColor3f(0.8,0.5,0);
    //glVertex3f(x[0],x[1],x[2]);
  }
  //glEnd();

  //glPopMatrix();


}

void draw_galactic_plane(const int num_grid)
{
  int i,k,ii;
  //COORD3D *ring1=new COORD3D [num_grid];
  COORD3D *ring2=new COORD3D [num_grid];
  COORD3D *ring3=new COORD3D [num_grid];

  //make_circle_galactic(3*1.02, 90, 360, ring1);
  make_circle_galactic(3*1.07, 90-1.5, num_grid, ring2);
  make_circle_galactic(3*1.07, 90+1.5, num_grid, ring3);
  
  double **mat_z=create_2d_grid<double>(3,3);  
  double **mat_ring=create_2d_grid<double>(3,1);
  double **mat_ring_new=create_2d_grid<double>(3,1);
  rot_z(90.0,mat_z);
  
  for(i=0;i<num_grid;i++)
  {
    for(k=0;k<3;k++) mat_ring[k][0]=ring2[i].r[k];
	AdotB(mat_z,mat_ring,mat_ring_new,3,3,1);
	for(k=0;k<3;k++) ring2[i].r[k]=mat_ring_new[k][0];
	
	for(k=0;k<3;k++) mat_ring[k][0]=ring3[i].r[k];
	AdotB(mat_z,mat_ring,mat_ring_new,3,3,1);
	for(k=0;k<3;k++) ring3[i].r[k]=mat_ring_new[k][0];
  }

  

  double p1[3],p2[3],p3[3],p4[3];
  double alpha=0.5;

  for(i=0;i<num_grid;i++)
  {
    for(k=0;k<3;k++) 
    {
      p1[k]=ring2[i].r[k];

      if(i==(num_grid-1)) ii=0;
      else ii=i+1;

      p2[k]=ring2[ii].r[k];
      p3[k]=ring3[ii].r[k];

      p4[k]=ring3[i].r[k];
    }


    glPushMatrix();
	glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.8,0.5,0,alpha);
    glBegin(GL_QUADS);
      glVertex3f(p1[0],p1[1],p1[2]);
      glVertex3f(p2[0],p2[1],p2[2]);
      glVertex3f(p3[0],p3[1],p3[2]);
      glVertex3f(p4[0],p4[1],p4[2]);
    glEnd();
    glPopMatrix();
  }

  //delete ring1;
  delete ring2;
  delete ring3;
  
  delete_2d_grid<double>(mat_ring,3,1);
  delete_2d_grid<double>(mat_ring_new,3,1);
  delete_2d_grid<double>(mat_z,3,3);
}



//tilt_plane [deg]
//rot_plane [deg]
void draw_sun_position(const double tilt_plane, const double rot_plane, const int num_grid)
{
  int i,j,k;

  COORD3D *ring=new COORD3D [num_grid];


  double **mat_z=create_2d_grid<double>(3,3);  
  double **mat_x=create_2d_grid<double>(3,3);  
  double **mat_zx=create_2d_grid<double>(3,3);
  double **mat_ring=create_2d_grid<double>(3,1);
  double **mat_ring_new=create_2d_grid<double>(3,1);

  rot_x(tilt_plane,mat_x);
  rot_z(rot_plane,mat_z);
  AdotB(mat_z,mat_x,mat_zx,3,3,3);

  make_circle_h(3*1.07, 90, num_grid, ring);

  for(i=0;i<num_grid;i++)
  {
    for(k=0;k<3;k++) mat_ring[k][0]=ring[i].r[k];
	AdotB(mat_zx,mat_ring,mat_ring_new,3,3,1);
	for(k=0;k<3;k++) ring[i].r[k]=mat_ring_new[k][0];
  }

  cout<<"Ring is prepared!"<<endl;

  double theta,phi;
  double x[3];
  double radius;

  radius=3*1.07;
  theta=90;
  phi=_CLOCK_SUN_*360.0/24.0;

  x[0]=radius*sin(deg2rad*theta)*cos(deg2rad*phi);
  x[1]=radius*sin(deg2rad*theta)*sin(deg2rad*phi);
  x[2]=radius*cos(deg2rad*theta);

  for(k=0;k<3;k++) mat_ring[k][0]=x[k];
  AdotB(mat_zx,mat_ring,mat_ring_new,3,3,1);
  for(k=0;k<3;k++) x[k]=mat_ring_new[k][0];



  int num_theta=100;
  int num_phi=200;

  COORD3D **sphere=new COORD3D * [num_phi];
  for(i=0;i<num_phi;i++) sphere[i]=new COORD3D [num_theta];

  draw_sky(0.2,num_phi,num_theta,sphere);

  cout<<"Sphere is prepared!"<<endl;


  for(i=0;i<num_phi;i++) for(j=0;j<num_theta;j++) for(k=0;k<3;k++) sphere[i][j].r[k]+=x[k];

  
  double p1[3],p2[3],p3[3],p4[3];
  int ii,jj;
  double alpha=0.3;

  glPushMatrix();
  glColor3f(0.5,0.8,0);


  glBegin(GL_LINES);
  for(i=0;i<num_grid;i++) glVertex3f(ring[i].r[0],ring[i].r[1],ring[i].r[2]);
  glEnd();

  glBegin(GL_LINES);
  glVertex3f(0,0,0);
  glVertex3f(x[0],x[1],x[2]);
  glEnd();

  for(i=0;i<num_phi;i++) for(j=0;j<num_theta;j++)
  {
    if(i==(num_phi-1)) ii=0;
    else ii=i+1;

    if(j==(num_theta-1)) jj=0;
    else jj=j+1;

    for(k=0;k<3;k++) p1[k]=sphere[i][j].r[k];
    for(k=0;k<3;k++) p2[k]=sphere[ii][j].r[k];
    for(k=0;k<3;k++) p3[k]=sphere[ii][jj].r[k];
    for(k=0;k<3;k++) p4[k]=sphere[i][jj].r[k];

	glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(1.0,0.0,0,alpha);
    glBegin(GL_QUADS);
      glVertex3f(p1[0],p1[1],p1[2]);
      glVertex3f(p2[0],p2[1],p2[2]);
      glVertex3f(p3[0],p3[1],p3[2]);
      glVertex3f(p4[0],p4[1],p4[2]);
    glEnd();
  }
  glPopMatrix();

  delete_2d_grid<double>(mat_ring,3,1);
  delete_2d_grid<double>(mat_ring_new,3,1);
  delete_2d_grid<double>(mat_z,3,3);
  delete_2d_grid<double>(mat_x,3,3);
  delete_2d_grid<double>(mat_zx,3,3);

  for(i=0;i<num_theta;i++) delete sphere[i];
  delete sphere;

  delete ring;

}


void draw_background1()
{
  int i,j;
  int num_theta=1000;
  int num_phi=2000;
  double dtheta=180.0/num_theta;
  double dphi=360.0/num_phi;
  double theta,phi;
  double x,y,z;
  hpint64 idpix;
  int nside=_nside_;

  double r=6*1.07;
  double rr,gg,bb;
  double ngal;
  
  int ncolor=256*256;

  reset_color();
  
  

  glPushMatrix();
  //double rot=_ID_STEP_*185.0/3600.0*360.0/24.0*1.6;
  //rotate_background();
  //glRotatef(rot, 0.0f, 1.0f, 0.0f);
  //glRotatef(-90,1.0f,0.0f,0.0f);

  Healpix_Map<double>  *pMap;
  pMap=&_map_GAL_;

  glBegin(GL_POINTS);
  for(i=0;i<pMap->Npix();i++)
  {
    ngal=(*pMap)[i];

    if(ngal<=-0.5) ngal=-0.5;
    if(ngal>=0.5) ngal=-0.5;
    ngal=(ngal+0.5)*ncolor;

    //color2rgb(ngal,rr,gg,bb);

    pix2ang_ring64(nside, i, &theta, &phi);

    double radius=0.2*deg2rad;
    pointing pt;
    vector<int> listpix;

    pt.theta=theta;
    pt.phi=phi;
/*
    pMap->query_disc_inclusive(pt,radius,listpix);
    double v=0.0;
    for(int ii=0;ii<listpix.size();ii++) v+=(*pMap)[listpix[ii]];
    v+=(*pMap)[i];
    v/=(listpix.size()+1.0);
    ngal=v;
*/
//    if(ngal<=-0.5) ngal=-0.5;
//    if(ngal>=0.5) ngal=-0.5;
    //ngal=(ngal+0.5)*ncolor;

    color2rgb(ngal,rr,gg,bb);
    glColor3f(rr,gg,bb);
/*
//if(listpix.size()>=3)
if((phi<5.0*PI/4.0)&&(phi>-PI/2.0))
{
    glBegin(GL_POLYGON);
    for(int ii=0;ii<listpix.size();ii++)
    {
      pix2ang_ring64(nside, listpix[ii], &theta, &phi);
      x=r*sin(theta)*cos(phi);
      y=r*sin(theta)*sin(phi);
      z=r*cos(theta);

      //glColor3f(rr,gg,bb);
      glVertex3f(x,y,z);
    }
    glEnd();
}
*/
//	if((phi<5.0*PI/4.0)&&(phi>-PI/2.0))
	//{

    x=r*sin(theta)*cos(phi);
    y=r*sin(theta)*sin(phi);
    z=r*cos(theta);
    glColor3f(rr,gg,bb);
    glVertex3f(x,y,z);
	//}

  }
  glEnd();

  glPopMatrix();
}

/*
void draw_bk()
{

  reset_color();

  double rot=0;//id_sphere*185.0/3600.0*360.0/24.0;

  glRotatef(rot, 0.0f, 1.0f, 0.0f);

  glRotatef(-90,1.0f,0.0f,0.0f);
  gluQuadricTexture(quad_bk,1);
  gluSphere(quad,6*1.07,100,100);

}



void draw_sun()
{
  glPushMatrix();
  glRotatef(-90,1.0f,0.0f,0.0f);
  glTranslatef(3.5,5.0,1.6);

  //gluQuadricTexture(quad,1);
  reset_color();

  GLUquadric *quad1 = gluNewQuadric();
  gluSphere(quad1,0.5,100,100);
  delete quad1;
  glPopMatrix();
}
*/


