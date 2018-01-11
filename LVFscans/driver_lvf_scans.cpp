#include "src/lvf_scans.h"


#include "src/global.h"
#include "src/systematics.h"


int main(int argc, char *argv[])
{
  _num_orbit_=1;//atof(argv[1]);
  
  last_mat=create_2d_grid<double>(3,3);
  last_pos=new double [3];
  last_camera_pointing=new double [3];

  //GLuint id_texture;

  TEXTURE_w_earth=1024;
  TEXTURE_h_earth=512;

  TEXTURE_w_sun=1024;
  TEXTURE_h_sun=512;

  TEXTURE_DATA_earth=load_jpg_X("/home/cfeng/computing/spherex/LVFscans/texture/earth.bmp",TEXTURE_w_earth, TEXTURE_h_earth);
  TEXTURE_DATA_sun=load_jpg_X("/home/cfeng/computing/spherex/LVFscans/texture/gal.bmp",TEXTURE_w_sun, TEXTURE_h_sun);
  //cout<<"Loaded textures"<<endl;
  //load_jpg1("/home/cfeng/computing/spherex/LVFscans/texture/gal.bmp",1024,512);

  Healpix_Map<double> planckmap143(1024,RING,SET_NSIDE);
  //load_fits_map_cxx("/home/cfeng/computing/spherex/LVFscans/output/HFI_SkyMap_143_2048_R1.10_nominal.fits",planckmap143);
  load_fits_map_cxx("/home/cfeng/computing/spherex/LVFscans/output/SPHEREx_sky.fits",_map_GAL_);
  //_map_GAL_.Import_degrade(planckmap143);

  //DEEP_pointings=load_pointings("/home/cfeng/computing/spherex/LVFscans/output/deep_pointings");
  DEEP_pointings=load_pointings("/home/cfeng/computing/spherex/LVFscans/output/deep_pointing_T5");
  //ROTATE_LVF_90(pos);
//cout<<DEEP_pointings[0][0]<<"\t"<<DEEP_pointings[0][1]<<endl;
//cout<<DEEP_pointings[6][0]<<"\t"<<DEEP_pointings[6][1]<<endl;

  dtr=init_detector_view(WIDTH_DETECTOR,HEIGHT_DETECTOR);

  prepare_systematics();

#ifdef _SHOW_DECTECTOR_HISTORY_
  DH=new DETECTOR_HISTORY;
  DH->num_dtr=0;
  DH->default_num_dtr=10000;

  DH->dtr_arr=new DETECTOR * [DH->default_num_dtr];
  for(int i=0;i<DH->default_num_dtr;i++) DH->dtr_arr[i]=init_detector_view(WIDTH_DETECTOR,HEIGHT_DETECTOR);
#endif



#ifdef _USE_OPENGL_
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);


  windW=800;//2560;
  windH=600;//1440;
  glutInitWindowSize(windW,windH);
  glutInitWindowPosition(0, 0);
  mw=glutCreateWindow("SphereX");
  glutSetWindow(mw);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


  glutCreateMenu(menu);
  glutAddMenuEntry("Show Background Galaxy", 1);
  glutAddMenuEntry("Don't Show Background Galaxy", 0);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  //quad_sun=init_texture(TEXTURE_w_sun, TEXTURE_h_sun, texture_sun, TEXTURE_DATA_sun);
  //glUseProgram(id_texture);
//  quad_earth=
init_texture(0,TEXTURE_w_earth, TEXTURE_h_earth, texture_earth, TEXTURE_DATA_earth);
//  quad_sun=
//init_texture(1,TEXTURE_w_sun, TEXTURE_h_sun, texture_sun, TEXTURE_DATA_sun);

  //cout<<"Quads for Textures"<<endl;

  init();
  glutDisplayFunc(get_position);
  glutTimerFunc(100, timer, 0);
  glutReshapeFunc(reshape);


  glutMouseFunc(mouseButton);
  glutMotionFunc(mouseMove);

  glutMainLoop();


#else
  //load_fits_map_cxx("/home/cmb/computing/spherex/LVFscans/output/wmap_band_imap_r9_9yr_W_v5.fits");
  //test_kid2parent();return 1;
/*
int len=5;
double arr[5];
arr[0]=2;
arr[1]=9;
arr[2]=3;
arr[3]=5.9;
arr[4]=1.9;

double a=get_list_minimum(len,arr);
cout<<a<<endl;throw;
*/
  init();
  get_original_scan();//get_grid(0,_map_GAL_, "grid_original");
  for(int i=0;i<1.5*365*8;i++)
  //for(int i=0;i<1.9e5;i++)
  {
    _ID_STEP_=i;
    //rotate_background();
    get_position(0);

    //_ID_STEP_=i;
  }
  //Analyze_NCP_flat();
  //OUTPUT_NCP_flat();
  //postprocessing();
#endif


  return 1;
}
