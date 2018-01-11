#include "lvf_scans.h"

#include "global.h"


void get_position(void)
{
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);

  glLoadIdentity();
  //gluLookAt(1,6,6,0,0,0,0,1,0);
  gluLookAt(CAMERA_X,CAMERA_Y,CAMERA_Z,0,0,0,0,1,0);
  //gluLookAt(0,0,7,0,0,0,0,1,0);
  //cout<<"Camera position: "<<CAMERA_X<<"\t"<<CAMERA_Y<<"\t"<<CAMERA_Z<<endl;


  glEnable(GL_DEPTH_TEST);

  glRotatef(-90,1,0,0);
  glRotatef(-90,0,0,1);

  glPushMatrix();
  earth_texture();
  glPopMatrix();

  //make_circle_galactic(3*1.02, 90, 360);
  //make_circle_galactic(3*1.02, 90-1.5, 360);
  //make_circle_galactic(3*1.02, 90+1.5, 360);

  if(FLAG_show_ncp) show_ncp(360);

  //draw_galactic_plane(360);

  //draw_sun_position(_SUN_angle_1_, _SUN_angle_2_, 360);


  if(_FLAG_show_axis_)
  {
    glPushMatrix();
    drawAxes(1);
    glPopMatrix();
  }
  else
  {
    glPushMatrix();
    drawAxes(1e-3);
    glPopMatrix();
  }

  if(_FLAG_show_galaxy_)
  {
    glPushMatrix();
    draw_background1();
    /*
    glEnable(GL_TEXTURE_2D);
    glBindTexture( GL_TEXTURE_2D, texture_bk);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    draw_bk();
    glDisable(GL_TEXTURE_2D);
    */
    glPopMatrix();
    //draw_sun();
  }

  glPushMatrix();
  get_position(0);
  glPopMatrix();

  glutSwapBuffers();
  
  
  //make_movie("scan/frame",windW,windH);
  
  _GLOBAL_CLOCK_+=1;
  //update_CAMERA();
  //update_cant_angle();
  //update_sat_theta();
  //update_sat_delta();
  _ID_STEP_++;
  _CLOCK_SUN_++;
  
  //cout<<"Main->"<<_ID_STEP_<<endl;


}
