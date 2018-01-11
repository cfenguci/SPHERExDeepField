#include "lvf_scans.h"

#include "global.h"

void show_ncp(const int num_grid)
{
  COORD3D *ring_small=new COORD3D [num_grid];
  COORD3D *ring_big=new COORD3D [num_grid];
  double alpha=0.6;
  int i,j,k,ii;
  double p1[3],p2[3],p3[3],p4[3];
  double pole[3];

  double far_radius=4*1.07;

  pole[0]=0;
  pole[1]=0;
  pole[2]=far_radius;

  make_circle_h(RADIUS_EARTH, POLE_BOUND, num_grid, ring_small);
  make_circle_h(far_radius, POLE_BOUND, num_grid, ring_big);

  glPushMatrix();
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0,1.0,1.0,alpha);
  
  for(i=0;i<num_grid;i++)
  {
    glBegin(GL_QUADS);

    if(i==(num_grid-1)) ii=0;
    else ii=i+1;

    for(k=0;k<3;k++)
    {
      p1[k]=ring_small[i].r[k];
      p2[k]=ring_small[ii].r[k];
      p3[k]=ring_big[ii].r[k];
      p4[k]=ring_big[i].r[k];
    }
    glVertex3f(p1[0],p1[1],p1[2]);
    glVertex3f(p2[0],p2[1],p2[2]);
    glVertex3f(p3[0],p3[1],p3[2]);
    glVertex3f(p4[0],p4[1],p4[2]);

    glEnd();
  }
  glPopMatrix();

  glPushMatrix();
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0,0.0,0.0,alpha);
  for(i=0;i<num_grid;i++)
  {
    glBegin(GL_TRIANGLES);

    if(i==(num_grid-1)) ii=0;
    else ii=i+1;

    for(k=0;k<3;k++)
    {
      p1[k]=ring_big[i].r[k];
      p2[k]=ring_big[ii].r[k];
      p3[k]=pole[k];
    }
    glVertex3f(p1[0],p1[1],p1[2]);
    glVertex3f(p2[0],p2[1],p2[2]);
    glVertex3f(p3[0],p3[1],p3[2]);
    glEnd();
  }
  
  glPopMatrix();

  delete ring_small;
  delete ring_big;
}
