#!/bin/bash
#probleme avec intel/17 sur un test
set -e

Untar(){
	echo "======== I untar my package ========"
#	cp $software-${version%%-*}.tar.gz $compile_dir/.
	cd $compile_dir
	tar -zxvf  /opt/software/tarballs/$software-${version%%-*}.tgz
}

Configure(){
	echo "======== I configure ========"
	cd $compile_dir/$software-${version%%-*}
	export CPATH=$MKLROOT/include/fftw:$MKLROOT/include:$CPATH
	./configure FC=mpif90 F90=mpif90 F77=mpif90 CC=mpicc CXX=mpicxx CXXFLAGS="-align -g -traceback -qoverride-limits -diag-disable:remark -O2 -mcmodel=medium -mkl=sequential -xCORE-AVX2 -fPIC -shared-intel -lfftw3xc_intel" CFLAGS="-align -g -traceback -qoverride-limits -diag-disable:remark -O2 -mcmodel=medium -mkl=sequential -xCORE-AVX2 -fPIC -shared-intel -lfftw3xc_intel" LDFLAGS="-align -g -traceback -qoverride-limits -diag-disable:remark -O2 -mcmodel=medium -mkl=sequential -xCORE-AVX2 -fPIC -shared-intel -lfftw3xc_intel" --prefix=$install_dir --disable-openmp --enable-fftw
}

Test(){
	echo "======== I test ========"
	export PATH=$install_dir/bin:$PATH
	export CPATH=$install_dir/include:$CPATH
	export FPATH=$install_dir/include:$FPATH
	export LIBRARY_PATH=$install_dir/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$install_dir/lib:$LD_LIBRARY_PATH
	export PLUMED_KERNEL=$install_dir/lib/libplumedKernel.so
#	make check
#	make test
}

Install(){
	echo "======== I install  ========"
	make -j 40
	make -j 40 install
}

Clean(){
	echo "======== I clean my directory ========"
	make clean
}

Deploy(){

    # if [ -d $install_dir ];then  
    # 	rm -rf $install_dir
    # fi
    # mkdir -p $install_dir
    
    # if [ -d $compile_dir ];then  
    # 	rm -rf $compile_dir
    # fi
    # mkdir -p $compile_dir

    Untar
    module list
    Configure
#    Test
    Install
#    Clean
}

Deploy
