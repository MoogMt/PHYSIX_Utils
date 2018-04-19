
COMPILATION
========================

Works with Gromacs 5.1.4 et Plumed 2.4.1

----------------------------------------
1- Compile your own FFTW3 LIBRARY:
----------------------------------------

-> Get archive from site, dearchive and go inside folder

http://www.fftw.org/download.html 

tar -xvf fftw3_archive

cd fftw3_folder

-> Preparing Install

fftw_folder=$(pwd)"/install_fftw" 

mkdir $fftw_folder 
-> Configure with float, mpi, shared version and in install folder

./configure --prefix=$folder --enable-shared --enable-sse2 --enable-float --enable-mpi

-> Compiling f

make -j 8

make install -j 8

-> Moving back

cd .. 

----------------------------------------

----------------------------------------
2 - Getting GROMACS
----------------------------------------

git clone https://github.com/gromacs/gromacs

----------------------------------------

--------------------------------------
3 - Optionnal - PLUMED PATCHING
--------------------------------------

-> Get archive from plumed repository

git clone https://github.com/plumed/plumed2

-> Moving into folder

cd plumed2

-> Configuring install

folder_plumed=$(pwd)"install_plumed"

./configure --prefix=$folder

-> Installing

make -j 8

make -j 8 install

-> Sourcing plumed-path

source sourcememe.sh

-> Moving back into gromacs folder

cd ../gromacs

-> Patching 

plumed-patch -p --shared 

--------------------------------------

--------------------------------------
3 - Bis - Installing GROMACS
--------------------------------------

-> Moving to gromacs folder

cd gromacs

-> Making compilation and installation folder

mkdir build
mkdir install_build
grmx_install=$(pwd)"/grmx_install"

-> configuring

fftw_include=$fftw_folder"/include"

fftw_lib=$fftw_folder"/lib"

cmake .. -DGMX_MPI=on -DCMAKE_INSTALL_PREFIX=$grmx_install -DGMX_FFT_LIBRARY=fftw3 -DCMAKE_PREFIX_PATH=$fftw_lib-DFFTWF_INCLUDE_DIR=$fftw_folder

-> Compiling

make -j 8

/!\ Everything explodes here...

make -j 8 install

--------------------------------------
