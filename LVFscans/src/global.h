#include "lvf_scans.h"
#include "raj_code.h"

extern GLint windW, windH;
extern double _num_orbit_;
extern int _ID_STEP_;
extern int _GLOBAL_CLOCK_;
extern hpint64 _nside_,_nside_256_,_idpix_;
extern Healpix_Map<double> _map_LVF_scan_;

extern int id_sphere;
extern int _cnt_days_;



extern GLUquadric *quad_earth;
extern GLUquadric *quad_sun;

extern GLuint texture_earth;
extern GLuint texture_sun;


extern int TEXTURE_w_earth,TEXTURE_h_earth;
extern unsigned char *TEXTURE_DATA_earth;

extern int TEXTURE_w_sun,TEXTURE_h_sun;
extern unsigned char *TEXTURE_DATA_sun;

extern GLUquadric *quad_bk;
extern GLuint texture_bk;
extern unsigned char *_image_bk_;

extern int _ID_FOV_record_;
extern Healpix_Map<double> _map_GAL_;
extern Healpix_Map<double> _map_GAL_256_;
extern double _max_ra_,_max_dec_;
extern double _min_ra_,_min_dec_;

extern Healpix_Map<double> _map_NAIVE_;
extern Healpix_Map<hpint64> _index_mapping_;

extern int _rand_seed_;
extern planck_rng _rng_;

extern GLuint mw,subw;



//record the camera view
extern DETECTOR *dtr;
extern DETECTOR_HISTORY *DH;
extern int CURRENT_STRIP_ID;
extern int NUM_STRIP;

extern int _FLAG_show_galaxy_;
extern int _FLAG_show_axis_;

extern double CAMERA_X,CAMERA_Y,CAMERA_Z;
extern double last_camera_x,last_camera_y,last_camera_z;
extern double DefaultCameraX, DefaultCameraY, DefaultCameraZ;

extern double **last_mat,*last_pos,*last_camera_pointing,last_camera_theta,last_camera_phi;

extern double INIT_THETA,INIT_PHI;
extern double TEST_THETA,TEST_NOD;
extern double TEST_CANT;

extern LVF_STATUS LVF_STATUS_kid;

extern double _SUN_angle_1_, _SUN_angle_2_;

extern int FLAG_show_ncp;

extern int _CLOCK_SUN_;


extern struct_lvf_track *LVF_track;

extern int Count_NCP_pointings;

extern int Count_NCP_frames;

#include "util_gridxy.h"
extern GRIDXY **GDXY;


extern double **DEEP_pointings;

extern struct_data_table *dt_dark_current,*dt_thermo,*dt_zodi_A,*dt_zodi_B;
extern double **mat_dark_current,**mat_thermo,**mat_zodi_A,**mat_zodi_B;
extern double **mat_dark_current_LOW,**mat_thermo_LOW,**mat_zodi_A_LOW,**mat_zodi_B_LOW;


extern GRIDXY *image_ncp_s,*image_ncp_n,*image_ncp_zodi,*image_ncp_dark,*image_ncp_thermo;
