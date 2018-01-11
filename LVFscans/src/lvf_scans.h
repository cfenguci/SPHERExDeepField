#ifndef lvf_scans_h
#define lvf_scans_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "math.h"

#include "string.h"

#define PI 3.1415926535
#define deg2rad (PI/180.0)
#define MATTYPE_ZERO 0.00


#define LIGHT_SPEED 299792458 //m/s
#define _G_ 6.67428e-11 //m^3*kg^{-1}*s^{-2}
#define M_earth 5.972e24 //kg
#define R_earth 6.4e6 //m


#define WIDTH_DETECTOR 7.0
#define HEIGHT_DETECTOR 3.5
#define PIX_REDUC 32
#define DISPLAY_REDUC 1

#define POLE_BOUND 10.0
#define POLE_BOUND_GRID 7.0
#define FLAT_PATCH 0.4
#define FLAT_SPACE 0.001//0.001

#define NUM_STEPPING 1
#define RADIUS_EARTH 0.5

//#define FAR_VIEW (10000*1.07)
#define FAR_VIEW (3*1.07)

#include "healpix_include.h"

#include <GL/glut.h>


//#define _rotate_lvf_

//#define _USE_OPENGL_
//#define _SHOW_DECTECTOR_HISTORY_
//#define _MOVIE_MODE_

struct COORD
{
  double ra,dec;
};

struct COORD3D
{
  double r[3];
};

struct DETECTOR
{
  int numx,numy;
  double pix,xmin,ymin;
  double **fov_image;
  COORD **fov_coord;

  double **fov_coord3d_x;
  double **fov_coord3d_y;
  double **fov_coord3d_z;

  double omega_theta,omega_phi;

  double strip_width;
  double delta_t;
};

struct DETECTOR_HISTORY
{
  DETECTOR **dtr_arr;	
  int num_dtr;
  int default_num_dtr;

};


struct LVF_STATUS
{
  int flag_show_axis;
  int flag_show_background;
  int step;
};

struct pair_mat
{
  double **mat_theta,**mat_phi;
};

struct struct_lvf_track
{
  int num_step;
  int num_grid;
  COORD3D **STEP_RING;
  int current_num_rings;
};
//instrument.cpp
struct_lvf_track * new_LVF_track(const int num_step, const int num_grid);
void delete_LVF_track(struct_lvf_track *p, const int num_step, const int num_grid);
void insert_LVF_track(COORD3D *, COORD3D *);



//coor_trans.cpp
void equ2gal(const double ra, const double dec, double &l, double &b);
void gal2equ(const double l, const double b, double &ra, double &dec);


//utils.cpp
template<typename MATTYPE> MATTYPE **create_2d_grid(const int nx, const int ny)
{
  int i,j;
  MATTYPE **mat=(MATTYPE **)calloc(nx,sizeof(MATTYPE *));
  for(i=0;i<nx;i++) mat[i]=(MATTYPE *)calloc(ny,sizeof(MATTYPE));
  //if(MATTYPE!=COORD)
  for(i=0;i<nx;i++) for(j=0;j<ny;j++) mat[i][j]=(MATTYPE) MATTYPE_ZERO;
  //else
  //{
  //  for(i=0;i<nx;i++) for(j=0;j<ny;j++) {mat[i][j].ra=0;mat[i][j].dec=0;}
  //}
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


void vecA_cross_vecB(double *A, double *B, double *C);
void AdotB(double **A,double **B, double **C, const int dimA,const int dim, const int dimB);
void AdotB(double **A, double multi, const int dimA, const int dimB);

double det33(double **A);
void rot_x(const double theta,double **);
void rot_y(const double theta,double **);
void rot_z(const double theta,double **);
int sphere_normal_vec(const double theta, const double phi, double *x);

void get_matW(double *u, double **matW);
void get_arbitrary_rot(const double phi, double *u, double **mat);
double get_vec_len(double *vec);
void mat_cp(double **, double **);
void vec_cp(double *,double *);

COORD **create_2d_grid_coord(const int nx, const int ny);
void delete_2d_grid_coord(COORD ** mat, const int nx, const int ny);
COORD ***create_3d_grid_coord(const int nx, const int ny, const int nz);
void insert_detector(DETECTOR *p);

//lvf_scans.cpp
void make_circle(const double radius, const double phi, const int num_grid,const double r, const double g, const double b);
void make_circle_h(const double radius, const double theta, const int num_grid,const double r, const double g, const double b);
void make_circle_h(const double radius, const double theta, const int num_grid,COORD3D *ring);
void make_circle_v(const double radius, const double phi, const int num_grid,COORD3D *ring);
void generate_pointing(const double, const double theta_0, const double phi_0,const double omega_theta, const double,const double,double *,double *, double &, double &);
void generate_pointing_slews(const double, const double theta_0, const double phi_0, const double, const double, const double, const int id_slews, const double, double *x,double &theta_out, double &phi_out);

void spherical_2_cartesian(const double x, const double y, const double z, double &theta, double &phi);

void reset_color();
void draw_sphere(GLUquadric *tmp_quad);
void earth_texture();
void draw_sky(const double,const int, const int,COORD3D **);


void get_LVF_image(const int idi,const int idj,const int idpix, const double dec, const double ra, DETECTOR *dt);

void show_sub_window(const int numx, const int numy, double *xx, double *yy, double **vv);
void FOV_image(const int, DETECTOR *dt,double **rot_LVF, double *xpos_LVF);

void init();
void get_position_sub(void);

void update_sat_theta();
void update_sat_delta();
void update_cant_angle();
void update_CAMERA();


void get_strip(const int id_strip, const int, double &ymin, double &ymax);
DETECTOR * init_detector_view(const double, const double);
void delete_detector_view();



void get_pole(const int);

//event_handleX.cpp
void get_position(void);
void get_position(const int);



void mouseMove(int xpos, int ypos);
void mouseButton(int button, int state, int xpos, int ypos) ;
void timer(int v);
void timer1(int v);
void reshape(GLint w, GLint h);
void resize_sub(int width, int height);

void make_movie(const int, string name_movie, const int , const int );
void make_movie(string name_movie, const int , const int );


void load_fits_map_cxx(string file);
void load_fits_map_cxx(string file, Healpix_Map<double> &HealpixDATA);

void write_grid_data(DETECTOR *dtr);

void postprocessing();
void postprocessing(string file);



//backgroud_galaxy.cpp
void test_color();
void color2rgb(const int c, double &rr, double &gg, double &bb);
void get_map_pix_stat();
void rotate_background();
hpint64 get_index_mapping_parent2kid(hpint64 id_parent);
hpint64 get_index_mapping_kid2parent(hpint64 id_kid);
double *RING_position2vec(const double l, const double b);
void make_circle_galactic(const double radius, const double theta, const int num_grid, COORD3D *);
void draw_galactic_plane(const int);
void draw_sun_position(const double tilt_plane, const double rot_plane, const int num_grid);
void draw_background1();


//axis.cpp
void Arrow(GLdouble x1,GLdouble y1,GLdouble z1,GLdouble x2,GLdouble y2,GLdouble z2,GLdouble D);
void drawAxes(GLdouble length);

//write_bmp.cpp
void WriteBMP(const char *fname, int w,int h,unsigned char *img);


//menu.cpp
void menu(int item);

//object_texture.cpp
unsigned char * load_jpg_X(char *filename, const int width, const int height);
void init_texture(const int, const int width, const int height, GLuint tt, unsigned char * TEXTURE_DATA);

//spacecraft.cpp
void show_spacecraft_z(const double, const double,double **,double *);
void get_normal_plane_axis(double **rot_LVF, double *xpos_LVF, double center_LVF[3],double axis_x[3], double axis_y[3]
,double center_LVF_3d[3],double axis_x_3d[3], double axis_y_3d[3]);

void get_FOV_center(double **axis_spacecraft_z, double &theta_camera_pointing, double &phi_camera_pointing);
void SHOW_normal_plane_axis(const double ratio, double center_LVF_3d[3],double axis_x_3d[3], double axis_y_3d[3]);
void get_normal_plane_trans_mat(const double, double center_LVF_3d[3],double axis_x_3d[3], double axis_y_3d[3],double **mat_normal);
void get_fov_plane_axis(const double len, double **mat_LVF_pointing,double *xpos_LVF, double **axis_spacecraft_x,double **axis_spacecraft_y,double **axis_spacecraft_z);
void SHOW_fov_plane_axis(double center_LVF_3d[3],double **axis_spacecraft_x,double **axis_spacecraft_y,double **axis_spacecraft_z);
void SHOW_vh_circles(const double strip_width, const int, const double last_camera_phi, double **axis_spacecraft_z);

void get_instrument(double xFOV[4],double yFOV[4],double zFOV[4], double *xpos_LVF, double **rot_LVF);
void SHOW_camera(const double r, double xFOV[4],double yFOV[4],double zFOV[4], double *xpos_LVF, double **mat_SKY_pointing,double **mat_slews);

void get_detector_pixel(const double r, double rFOV[3], double *xpos_LVF, double **mat_SKY_pointing,double **mat_slews, double rFOV_sky[3]);
void record_detector_pixel(const double r, double *xpos_LVF, double **mat_SKY_pointing,double **mat_slews);

pair_mat * instrument_maneuver();
void delete_maneuver_mat(pair_mat *pm);

//spacecraft_patch.cpp
void show_spacecraft_z_patch(const double, const double,double **,double *);

//powspec.cpp
void make_mask(const double theta_0, Healpix_Map<double> &mask);
void get_naive_map_powspec();


//instrument.cpp
void copy_vec2vec(double *vec_object, double *vec_object_new);
void copy_vec2mat(double *vec_object, double **mat_object);
void copy_mat2vec(double **mat_object, double *vec_object);
void move_object(double **mat, double *pos, double *object, double *object_new);
void instrument(double *center_LVF, double **mat_SKY_pointing,pair_mat *);

double Shoot_LVF(double **mat_inst, double **mat_boresight, double *, double height_LVF, double *);

int get_LVF_band(const double x, const double y);

void shoot_deep(const int,double **,const double theta, const double phi);
void Shoot_LVF_fast(double **mat_inst, double **mat_boresight, double *, double height_LVF, double *);

void SHOW_detector2boresight(double *pos);
void get_LVF_center_circles(double *pos);

void ROTATE_LVF_90(double *pos);


//show_ncp.cpp
void show_ncp(const int num_grid);


#include "flat_analysis.h"
#include "flat_ncp.h"

#include "load_pointings.h"

#endif
