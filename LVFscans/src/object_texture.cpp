#include "lvf_scans.h"

#include "global.h"





unsigned char * load_jpg_X(char *filename, const int width, const int height)
{
  FILE * file;
  file = fopen( filename, "rb" );
  if (file == NULL) return 0;
  unsigned char * texture_data = (unsigned char *)malloc(width * height * 3 *sizeof(unsigned char));
  fread(texture_data, width * height * 3, 1, file);
  fclose(file);

  return texture_data;
}

void init_texture(const int which,const int width, const int height, GLuint tt, unsigned char * TEXTURE_DATA)
{
  glGenTextures(1, &tt);

  glBindTexture(GL_TEXTURE_2D, tt);

  glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_SRC_COLOR);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, TEXTURE_DATA);

  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );

  gluBuild2DMipmaps(GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, TEXTURE_DATA);


if(which==0)
{
  quad_earth = gluNewQuadric();
  glEnable(GL_TEXTURE_2D);
  glBindTexture( GL_TEXTURE_2D, tt);
  gluQuadricTexture(quad_earth,1);
}
else
{

  quad_sun = gluNewQuadric();
  //glEnable(GL_TEXTURE_2D);
  //glBindTexture( GL_TEXTURE_2D, tt);
  //gluQuadricTexture(quad_sun,1);
}



}

