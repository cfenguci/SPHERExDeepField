#include "lvf_scans.h"

#define RADPERDEG 0.0174533

void Arrow(GLdouble x1,GLdouble y1,GLdouble z1,GLdouble x2,GLdouble y2,GLdouble z2,GLdouble D)
{
  double x=x2-x1;
  double y=y2-y1;
  double z=z2-z1;
  double L=sqrt(x*x+y*y+z*z);

    GLUquadricObj *quadObj;

    glPushMatrix ();

      glTranslated(x1,y1,z1);

      if((x!=0.)||(y!=0.)) {
        glRotated(atan2(y,x)/RADPERDEG,0.,0.,1.);
        glRotated(atan2(sqrt(x*x+y*y),z)/RADPERDEG,0.,1.,0.);
      } else if (z<0){
        glRotated(180,1.,0.,0.);
      }

      glTranslatef(0,0,L-4*D);

      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder(quadObj, 2*D, 0.0, 4*D, 32, 1);
      gluDeleteQuadric(quadObj);

      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluDisk(quadObj, 0.0, 2*D, 32, 1);
      gluDeleteQuadric(quadObj);

      glTranslatef(0,0,-L+4*D);

      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder(quadObj, D, D, L-4*D, 32, 1);
      gluDeleteQuadric(quadObj);

      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluDisk(quadObj, 0.0, D, 32, 1);
      gluDeleteQuadric(quadObj);

    glPopMatrix ();

}
void drawAxes(GLdouble length)
{
  GLdouble D;
  D=0.01;

  reset_color();

  glPushMatrix();
  glRotatef(90,1,0,0);
  glRotatef(-90,0,0,1);
  glRotatef(180,0,0,1);
  glColor3f(0.5,0,0);
  glTranslatef(0,0,0);
  Arrow(0,0,0,length,0,0, D);
  glPopMatrix();

  glPushMatrix();
  glRotatef(90,1,0,0);
  glRotatef(-90,0,0,1);
  glRotatef(180,0,0,1);
  glColor3f(0,0.5,0);
  glTranslatef(0,0,0);
  Arrow(0,0,0,0,-length,0, D);
  glPopMatrix();

  glPushMatrix();
  glRotatef(90,1,0,0);
  glRotatef(-90,0,0,1);
  glColor3f(0,0,0.5);
  glTranslatef(0,0,0);
  Arrow(0,0,0, 0,0,-length, D);
  glPopMatrix();

  reset_color();
}
