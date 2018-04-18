=== COMPILATION STANDARD ===

----------------------------------------
1- Compile your own FFTW3 LIBRARY:
----------------------------------------

-> Get archive from site, dearchive and go inside folder
http://www.fftw.org/download.html 
tar -xvf fftw3_archive
cd fftw3_folder

-> Configure Install
fftw_folder=$(pwd)"/install" # install folder
mkdir $fftw_folder           # making folder
# Configure with float, mpi, shared version and in install folder
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
/!\ Can't crack it yet /!\ 
/!\ Something is wrong somewhere in plumed...
-> Get archive from plumed repository
git clone https://github.com/plumed/plumed2
-> Moving into folder
cd plumed2
-> Configuring install
folder_plumed=$(pwd)"install"
./configure --prefix=$folder
-> Installing
make -j 8
make -j 8 install
-> Moving back
cd ..
--------------------------------------

--------------------------------------
3 - Bis - Installing GROMACS
--------------------------------------
--------------------------------------
