#include "lvf_scans.h"
#include "raj_code.h"

GLint windW, windH;
double _num_orbit_;
int _ID_STEP_=0;
int _GLOBAL_CLOCK_=0;
hpint64 _nside_=1024,_nside_256_=1024,_idpix_;
Healpix_Map<double> _map_LVF_scan_(_nside_,RING,SET_NSIDE);

int id_sphere=0;
int _cnt_days_=0;
//unsigned char *_IMAGE_EARTH_;
//int _IMAGE_WIDTH_;
//int _IMAGE_HEIGHT_;



GLUquadric *quad_earth;
GLUquadric *quad_sun;

GLuint texture_earth;
GLuint texture_sun;


int TEXTURE_w_earth,TEXTURE_h_earth;
unsigned char *TEXTURE_DATA_earth;

int TEXTURE_w_sun,TEXTURE_h_sun;
unsigned char *TEXTURE_DATA_sun;



GLUquadric *quad_bk;
GLuint texture_bk;
unsigned char *_image_bk_;

int _ID_FOV_record_=0;
Healpix_Map<double> _map_GAL_(_nside_,RING,SET_NSIDE);
Healpix_Map<double> _map_GAL_256_(_nside_256_,RING,SET_NSIDE);
double _max_ra_=-1e100,_max_dec_=-1e100;
double _min_ra_=1e100,_min_dec_=1e100;

Healpix_Map<double> _map_NAIVE_(_nside_,RING,SET_NSIDE);
Healpix_Map<hpint64> _index_mapping_(_nside_,RING,SET_NSIDE);

int _rand_seed_=0;
planck_rng _rng_(_rand_seed_);

GLuint mw,subw;



//record the camera view
DETECTOR *dtr;
DETECTOR_HISTORY *DH;
int CURRENT_STRIP_ID=0;
int NUM_STRIP=48;

int _FLAG_show_galaxy_=0;
int _FLAG_show_axis_=1;

double CAMERA_X,CAMERA_Y,CAMERA_Z;
double last_camera_x,last_camera_y,last_camera_z;
//double DefaultCameraX=0.6, DefaultCameraY=9, DefaultCameraZ=0;
double DefaultCameraX=0.6, DefaultCameraY=7, DefaultCameraZ=0;

double **last_mat,*last_pos,*last_camera_pointing,last_camera_theta,last_camera_phi;

//(-8,20,8)
double INIT_THETA=0,INIT_PHI=0;
double TEST_THETA=0,TEST_NOD=0;
double TEST_CANT=0.0;


LVF_STATUS LVF_STATUS_kid;

double _SUN_angle_1_=100, _SUN_angle_2_=271;

int _CLOCK_SUN_=0;

int FLAG_show_ncp=0;


struct_lvf_track *LVF_track;

int Count_NCP_pointings;

int Count_NCP_frames;

#include "util_gridxy.h"
GRIDXY **GDXY;


double **DEEP_pointings;


struct_data_table *dt_dark_current,*dt_thermo,*dt_zodi_A,*dt_zodi_B;
double **mat_dark_current,**mat_thermo,**mat_zodi_A,**mat_zodi_B;
double **mat_dark_current_LOW,**mat_thermo_LOW,**mat_zodi_A_LOW,**mat_zodi_B_LOW;
GRIDXY *image_ncp_s,*image_ncp_n,*image_ncp_zodi,*image_ncp_dark,*image_ncp_thermo;
