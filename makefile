CC=g++
CFLAGS=-Wall -I.
DEPS = utils.h atom.h cell.h cutoff.h contact_matrix.h xyz.h  histogram.h sim.h pdb.h 
OBJ = read.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

co2_analysis: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -rf *.o 
