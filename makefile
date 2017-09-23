CC=g++
CFLAGS=-I.
DEPS = contact_matrix.h cell.h atom.h xyz.h utils.h
OBJ = read.o contact_matrix.o cell.o atom.o xyz.o utils.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

read: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
