# -DMPI       : turns on MPI part of the codes
# -DFOURBYTES : increases precision of integer representation of piv from 2 to 4 bytes
#               (it can be very slow without sizably improving the results!)

#fdflags =  -fbounds-check #debug
fcflags =   -c -O3 -cpp -DMPI #compile
flflags = -O3 #link
f90 = mpif90 -ffree-line-length-none -frecord-marker=4 -fno-range-check $(fdflags)
#--------------------------------------------------------------------------
build = ./build
src   = ./src
exec    = piv_clustering.x
#--------------------------------------------------------------------------
own     = get_lines read_blocks real2shortint shortint2real sprint print_network distance counting_sort main 
mod     = clustering_mod
lapack  = dlasrt2 lsame xerbla

objects = $(mod:%=$(build)/%.o) $(own:%=$(build)/%.o) $(lapack:%=$(build)/%.o) 

#--------------------------------------------------------------------------
#EXECUTABLE
#--------------------------------------------------------------------------
all: $(exec)

$(exec): $(objects)
	$(f90) $(flflags) -o $@ $^

#--------------------------------------------------------------------------
#OBJECTS (when any .o has to be made; this is how we do it)
#--------------------------------------------------------------------------
$(build)/%.o:$(src)/%.f90
	$(f90) $(fcflags) -o $@ $< 

#--------------------------------------------------------------------------
#MISC
#--------------------------------------------------------------------------
clean:
	rm -f  $(objects) $(mod:%=%.mod)

cclean:
	rm -f  $(objects) $(exec) $(mod:%=%.mod)
