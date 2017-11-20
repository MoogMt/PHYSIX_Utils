CPP=g++
CFLAGS=-std=gnu++11 -O3 -W -Wall -ansi -pedantic -I.

DEPS = utils.h atom.h cell.h cutoff.h contact_matrix.h xyz.h  histogram.h sim.h pdb.h molecules.h lut.h bonds.h

OBJP = pressure.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJNN = NN.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJmols = make_molecules.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJangles = angles.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJhistratio = ratio_counting.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o 						

%.o: %.c $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)

pressure: $(OBJP)
	$(CPP) -o $@ $^ $(CFLAGS)
NN: $(OBJNN)
	$(CPP) -o $@ $^ $(CFLAGS)
molecule: $(OBJmols)
	$(CPP) -o $@ $^ $(CFLAGS)	
angles: $(OBJangles)
	$(CPP) -o $@ $^ $(CFLAGS)
ratio: $(OBJhistratio)
	$(CPP) -o $@ $^ $(CFLAGS)	

clean:
	rm -rf *.o 



