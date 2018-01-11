WD=cfeng
HP_version=3.11


HPDIR=/home/$(WD)/computing/Healpix_$(HP_version)

INC=-I$(GSLDIR)/include $(INC_MOD)
LIB=-L$(GSLDIR)/lib -lgsl -lgslcblas $(LIB_GIF) \
    -L/home/$(WD)/computing/lapack-3.5.0 -llapack  -L/home/$(WD)/computing/lapack-3.5.0/BLAS/SRC -lblas_LINUX $(LIB_CAMB)


INC_CFITS=-I$(UCDIR)/CFITS/include
LIB_CFITS=-L$(UCDIR)/CFITS/lib -lcfitsio -L$(GFORTRAN) -lgfortran $(LIB)

INC_WCS=-I$(WCS)/include/wcslib
LIB_WCS=-L/home/cmb/computing/correlations/gamma_data/Fermi_LAT_app/wcslib-4.25.1/C -lwcs-4.25.1


DIR_cfits=$(CFITS)
DIR_healpix=$(HPDIR)#/home/cmb/computing/Healpix_3.20
INC_healpix=-I$(DIR_healpix)/include -I$(DIR_healpix)/src/f90/mod -I$(DIR_healpix)/src/cxx/generic_gcc/include
LIB_healpix= -L$(DIR_healpix)/lib -lgif -lhealpix -lhpxgif -lchealpix -lsharp_healpix_f \
	     -L$(DIR_healpix)/src/cxx/generic_gcc/lib -lhealpix_cxx  -lcxxsupport -lsharp -lfftpack -lc_utils \
	     -L$(DIR_cfits)/lib -lcfitsio

SRC_HP=/home/$(WD)/computing/Healpix_$(HP_version)/src/cxx/Healpix_cxx/healpix_base.cc#/home/cmb/computing/Healpix_3.20/src/cxx/build.generic_gcc/Healpix_cxx/healpix_base.cc


#########################################
INC_JPEG=-I$(JPEG)/include
LIB_JPEG=-L$(JPEG)/lib -ljpeg

VISUAL=/home/cfeng/computing/spherex/LVFscans/jpeg

SRC_JPEG=$(VISUAL)/example.c
OBJ_JPEG=$(VISUAL)/example.o
#########################################



GCC=gcc -std=gnu99
CFLAGS_GCC=$(INC_WCS) $(INC_CFITS)
LDFLAGS_GCC=$(LIB_WCS) $(LIB_CFITS) -lm

CC=g++
CFLAGS=-c -Wall $(INC_healpix) $(INC_CFITS) 
LDFLAGS= $(LIB_healpix) $(LIB_CFITS) -lm -lGL -lGLU -lglut


SOURCES= src/lvf_scans.cpp src/coor_trans.cpp src/utils.cpp src/backgroud_galaxy.cpp src/axis.cpp src/write_bmp.cpp src/menu.cpp src/object_texture.cpp \
src/spacecraft.cpp src/spacecraft_patch.cpp src/powspec.cpp \
src/event_handle.cpp src/event_handle_instance.cpp src/global.cpp \
src/instrument.cpp src/show_ncp.cpp \
src/raj_code.cpp src/flat_analysis.cpp \
src/util_gridxy.cpp src/flat_ncp.cpp \
src/load_pointings.cpp src/generate_pointings.cpp \
src/systematics.cpp src/bindata.cpp

SOURCES_WCS=src/get_fits_image.c src/get_fits_image_3d.c src/get_header.c src/wcs_pixel2radec.c src/wcs_coord.c 

OBJECTS=$(SOURCES:.cpp=.o) $(OBJ_HP)
#OBJECTS_WCS=$(SOURCES_WCS:.c=.o)

EXECUTABLE=bin/lvfscans bin/rings bin/postprocessing bin/raj bin/generate_pts bin/generate_ncp_noise bin/flat_image

all: $(OBJECTS) $(EXECUTABLE) 

bin/lvfscans: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_lvf_scans.cpp $(LDFLAGS) -o $@

bin/rings: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_rings.cpp $(LDFLAGS) -o $@

bin/postprocessing: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_postprocessing.cpp $(LDFLAGS) -o $@

bin/raj: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_raj.cpp $(LDFLAGS) -o $@

bin/generate_pts: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_generate_pointings.cpp $(LDFLAGS) -o $@

bin/generate_ncp_noise: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_generate_NCP_noise.cpp $(LDFLAGS) -o $@


bin/flat_image: $(OBJECTS)
	$(CC) $(INC_healpix) $(INC_CFITS) $(SRC_HP)  $(OBJECTS) driver_flat_image.cpp $(LDFLAGS) -o $@


#bin/use_fits: $(OBJECTS_WCS)
#	$(GCC) $(CFLAGS_GCC) $(OBJECTS_WCS) tfitshdr.c -o $@ $(LDFLAGS_GCC)

#bin/fermilat: $(OBJECTS_WCS)
#	$(GCC) $(CFLAGS_GCC) $(OBJECTS_WCS) driver_fermi_lat.c -o $@ $(LDFLAGS_GCC)

#bin/jpeg: $(OBJ_JPEG)
#	$(CC) $(FLAGS) $(INC) $(INC_JPEG) $(OBJ_JPEG) $(VISUAL)/jpeg.cpp -o bin/jpeg $(LIB_JPEG)




#$(OBJ_JPEG): %.o: %.c $(VISUAL)/*.h
#	cc  -c $(INC) $(FLAGS)  $(INC_JPEG) $(FLAGS) $< -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

#.c.o:
#	$(GCC) -c $(CFLAGS_GCC) $< -o $@

clean:
	rm src/*.o
	rm bin/*
