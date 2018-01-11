import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from os import listdir
from os.path import isfile, join
import subprocess
from scipy.interpolate import griddata
from math import pi  
import sys
#fig,ax=plt.subplots()
#cbar = fig.colorbar(cax)

whichgrid=sys.argv[1]
whichone=int(sys.argv[2])
whichtype=sys.argv[3]
vvmin=float(sys.argv[4])
vvmax=float(sys.argv[5])

nrows, ncols = 801, 801
FLAT_PATCH=0.4
FLAT_SPACE=0.001
dx=FLAT_SPACE
dy=FLAT_SPACE

fieldsize=8

nx=nrows/2-(fieldsize*pi/180.0)/dx
ny=ncols/2-(fieldsize*pi/180.0)/dy

def plot_sysmap_nproj():

  i=whichone
  fig,ax=plt.subplots()
 
  filename='framedata/daily/Flat_image_%s-%03d'%(whichgrid,i) 

  x,y,v,v1 = np.loadtxt(filename).T

  outputname='Flat_image_%s-%03d.png'%(whichgrid,i)

  grid = v.reshape((nrows, ncols))
  vx=(x-nrows/2.0)*dx*180/pi
  vy=(y-nrows/2.0)*dy*180/pi
  matx=vx.reshape((nrows, ncols))
  maty=vy.reshape((nrows, ncols))

  smatx=matx[nx:nrows-nx,ny:ncols-ny]
  smaty=maty[nx:nrows-nx,ny:ncols-ny]

  make_taper(smatx,smaty)

  subgrid=grid[nx:nrows-nx,ny:ncols-ny]
  submgrid = np.ma.masked_where(subgrid == 0, subgrid)  
  cax=ax.imshow(submgrid, extent=(smatx.min(), smatx.max(), smaty.max(), smaty.min()),interpolation='nearest', cmap=cm.get_cmap('rainbow'),vmin=vvmin,vmax=vvmax)
  cbar = fig.colorbar(cax)

  fig.savefig('framedata/daily/%s'%outputname)

def plot_sysmaps():

  i=whichone
  fig,ax=plt.subplots()
 
  filename='framedata/daily/Flat_image_%s-%03d'%(whichgrid,i);
  filename1='framedata/daily/Flat_image_%s-%03d-radial-template'%(whichgrid,i);
  filename2='framedata/daily/Flat_image_%s-%03d-radial-taper'%(whichgrid,i);
 

  x,y,v,v1 = np.loadtxt(filename).T


  xZ,yZ,vZ,v1Z = np.loadtxt(filename1).T;
  xZZ,yZZ,vZZ,v1ZZ = np.loadtxt(filename2).T;



  grid = np.zeros((nrows, ncols))
  outputname=''

  if whichtype=="mmp":
    v=vZZ*(v-vZ)
    outputname='Flat_image_%s-%03d_maskedmap.png'%(whichgrid,i)
    grid = v.reshape((nrows, ncols)).copy()

  if whichtype=="mp":
    outputname='Flat_image_%s-%03d.png'%(whichgrid,i)
    grid = v.reshape((nrows, ncols)).copy()


  if whichtype=="t":
    outputname='Flat_image_%s-%03d_template.png'%(whichgrid,i)
    grid = vZ.reshape((nrows, ncols)).copy()

  if whichtype=="m":
    outputname='Flat_image_%s-%03d_mask.png'%(whichgrid,i)
    grid = vZZ.reshape((nrows, ncols)).copy()


  


  vx=(x-nrows/2.0)*dx*180/pi
  vy=(y-nrows/2.0)*dy*180/pi
  matx=vx.reshape((nrows, ncols))
  maty=vy.reshape((nrows, ncols))

  smatx=matx[nx:nrows-nx,ny:ncols-ny]
  smaty=maty[nx:nrows-nx,ny:ncols-ny]

  subgrid=grid[nx:nrows-nx,ny:ncols-ny]


  submgrid = np.ma.masked_where(subgrid == 0, subgrid)

  
  cax=ax.imshow(submgrid, extent=(smatx.min(), smatx.max(), smaty.max(), smaty.min()),interpolation='nearest', cmap=cm.get_cmap('rainbow'),vmin=vvmin,vmax=vvmax)
  cbar = fig.colorbar(cax)

  fig.savefig('framedata/daily/%s'%outputname)
  

def make_taper(smatx,smaty):

  theta_1=4.0
  theta_2=5.5

  radius=np.sqrt(smatx**2+smaty**2)
  taper=np.zeros(np.shape(radius))
  taper[radius<=theta_1]=1.0
  taper[radius>theta_2]=0.0
  taper[(radius>=theta_1) & (radius<theta_2)]=cos(PI/2.0*(radius-theta_1)/(thetha_2-theta_1))

  fig,ax=plt.subplots()
  cax=ax.imshow(taper, extent=(smatx.min(), smatx.max(), smaty.max(), smaty.min()),interpolation='nearest', cmap=cm.get_cmap('rainbow'))
  cbar = fig.colorbar(cax)

  fig.savefig('deepfield_taper.pdf')
