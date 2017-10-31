CPP=g++
CFLAGS=-std=gnu++11 -O3 -W -Wall -ansi -pedantic -I.
DEPS = utils.h atom.h cell.h cutoff.h contact_matrix.h xyz.h  histogram.h sim.h pdb.h molecules.h lut.h bonds.h
OBJ = read.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o


%.o: %.c $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)

co2_analysis: $(OBJ)
	$(CPP) -o $@ $^ $(CFLAGS)

clean:
	rm -rf *.o 
