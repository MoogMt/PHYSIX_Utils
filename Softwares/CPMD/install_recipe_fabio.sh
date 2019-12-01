module purge
module load intel intelmpi
./scripts/configure.sh -DEST=exe LINUX-X86_64-INTEL-MPI-occigen

#   diff LINUX-X86_64-INTEL-MPI-occigen LINUX-X86_64-INTEL-MPI
#   22c22
#   <        FFLAGS='-O2 -cpp'
#   ---
#   >        FFLAGS='-O2'

cd exe
make -j 12

########## plumed:

export plumedir=/home/pietrucc/plumed-1.3-qe-pathcoordtable-cpmd4.1

cp $plumedir/common_files/patches/plumedpatch_cpmd-4.1.0.sh  exe

modify the patch to use a viable C++ compiler (on occigen: mpiicpc)

./plumedpatch_cpmd-4.1.0.sh -patch

make -j 12

