CC=g++
CFLAGS=-I.
DEPS = contact_matrix.h cell.h atom.h xyz.h utils.h histogram.h
OBJ = read.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

co2_analysis: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -rf *.o 
