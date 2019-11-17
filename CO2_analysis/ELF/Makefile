#-----------------------------------------------------------------------------#
# Makefile for cpmd2cube.x
# Configuration: Linux-PC-IFORT-FFTW
# Creation of Makefile: Jun 15 2018
# on Linux irene191 3.10.0-693.21.1.el7.x86_64 #1 SMP Fri Feb 23 18:54:16 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux
# Author: moogmath
#----------------------------------------------------------------------------
#
SHELL = /bin/sh
#
#--------------- Configuration
FC	= ifort
FFLAGS	=  -FR -D__IFC -DFFT_FFTW -O -tpp6 -w90 -w95 -pc64
LFLAGS	= -i-static -O -tpp6 -w90 -w95 -pc64
LIBS	= -I/ccc/cont003/home/gch0304/moogmath/fftw-2.1.5/install/include /ccc/cont003/home/gch0304/moogmath/fftw-2.1.5/install/lib/libfftw.a -L/ccc/cont003/home/gch0304/moogmath/fftw-2.1.5/install/lib
# uncomment and set only if needed
CC	= icc
CFLAGS	= -O2
#--------------- End of Configuration
##############

EXE	= cpmd2cube.x

OBJECTS	= \
    kinds.o periodic_table.o util.o grids.o cpmd2cube.o \
    grid_types.o atom_types.o readwrites.o \
    fftsg_lib.o fft_tools.o fftsg.o

CFGDEST = .
CFGMACH = Linux-PC-IFORT-FFTW
SRC=.

#############################

all:  $(EXE) trimcube.x cutcube.x

$(EXE): main.o $(OBJECTS)
	 $(FC) -o $(EXE) $(LFLAGS) main.o $(OBJECTS) $(LIBS)

cutcube.x: cutcube.o $(OBJECTS)
	 $(FC) -o $@ $(LFLAGS) cutcube.o $(OBJECTS) $(LIBS)

trimcube.x: trimcube.o
	$(CC) $(CFLAGS) -o $@ trimcube.o 
#############################
clean:
	rm -f *.o *.mod *.L *~

realclean: clean
	rm -f *.x

# Create a new makefile with new dependencies but retain the configuration.
reconfig:
	( cd $(SRC) ; $(SHELL) ./Configure -r -m -SRC=$(SRC) -DEST=$(CFGDEST) $(CFGMACH) ; cd $(CFGDEST) )

#############################
atom_types.o:  $(SRC)/atom_types.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

cpmd2cube.o:   $(SRC)/cpmd2cube.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

cutcube.o:   $(SRC)/cutcube.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

fft_tools.o:   $(SRC)/fft_tools.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

grid_types.o:  $(SRC)/grid_types.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

grids.o:       $(SRC)/grids.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

main.o:        $(SRC)/main.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

periodic_table.o: $(SRC)/periodic_table.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

readwrites.o:  $(SRC)/readwrites.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

util.o:       $(SRC)/util.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

kinds.o:      $(SRC)/kinds.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

fftsg.o:      $(SRC)/fftsg.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

fftsg_lib.o:  $(SRC)/fftsg_lib.F
	$(FC) -c $(FFLAGS) $(SRC)/$*.F

trimcube.o:	$(SRC)/trimcube.c
	$(CC) -c $(CFLAGS) $(SRC)/$*.c

# Dependencies
atom_types.o:  kinds.o \
		periodic_table.o

cpmd2cube.o:   atom_types.o \
		fft_tools.o \
		grid_types.o \
		grids.o \
		kinds.o \
		periodic_table.o \
		readwrites.o \
		util.o

cutcube.o:	kinds.o \
		periodic_table.o

fft_tools.o:   grid_types.o \
		fftsg_lib.o \
		kinds.o

grid_types.o:  kinds.o

grids.o:       grid_types.o \
		kinds.o \
		util.o

main.o:        cpmd2cube.o

periodic_table.o: kinds.o

readwrites.o:  atom_types.o \
		grid_types.o \
		kinds.o

util.o:       kinds.o

