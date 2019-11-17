#!/bin/bash
export target_libint=$1
if [[ $target_libint = "" ]]; then
	export target_libint=1.1.4
fi

source ./env
# sauvegarde du repertoire de lancement de compile.sh
set -e
export my_sweet_home=`pwd`

# La FFT est desormais ( a partir de la version 4.1 celle de la MKL


mkdir -p $compile_dir
# error if compile_dir path is too long : make[2]: execvp: /bin/sh: Argument list too long
echo "install_dir ...............",$install_dir
echo "compile_dir ...............",$compile_dir

Untar(){
  echo "************************************************************"
  echo "* Untar "
  echo "************************************************************"
	cd $compile_dir
  	#ls $HOME/tarballs/cp2k_$version.tar.gz
	pwd
	echo " ................................................................................................................................"
	echo " ......... Deploiement du tarball de CP2K version $version dans $compile_dir"
	echo " ................................................................................................................................"
	sleep 2
	rm -rf * # nettoyage de $compile_dir avant deploiement du tarball
  	tar -zxf /opt/software/tarballs/cp2k_$version.tar.gz
	find . -type d -name cp2k*4.1 
	find . -type d -name cp2k*4.1 | awk -f'/' {'print $2'}
	rep_asis=`find . -type d -name cp2k*4.1 | awk -f'/' {'print $2'}`
        repdir=`pwd`/$rep_asis
	echo "repdir =" $repdir
	ln -s $repdir/cp2k  $compile_dir/$version
	echo " affichage du contenu repertoire compile_dir/version"
	ls $compile_dir/$version
	cd -
}

Configure_and_prepare(){
  echo "************************************************************"
  echo "* Libraries internes à CP2K"
  echo "************************************************************"
# ..................................................................................................
# Configuration et compilation de libint 1.1.4
# librarie compilee en externe dependante de la version de intel ; nompi
# ..................................................................................................
# ..................................................................................................
# Configuration et compilation de  libxc 2.1.0
# 	en prodction pour l'instant cette lib XC n'est pas utilisee dans cp2k ( années 2015 et 2016)
# ..................................................................................................
# Configuration et compilation de tools/build_libsmm
# 	en prodction pour l'instant cette lib n'est pas utilisee dans cp2k ( années 2015 et 2016)
# ..................................................................................................

  if [[ $version == *"2.5.1"* ]] ; then
# Compilation du wrapper dans le sous repertoire tools/hfx_tools/libint_tools
	cd $compile_dir/$software-$version/tools/hfx_tools/libint_tools
	pwd
	ls
	icc -c libint_cpp_wrapper.o libint_cpp_wrapper.cpp -I$HOME/$machine/libraries/libint/$target_libint/$compiler/$compiler_version/nompi/include/libint -I$HOME/$machine/libraries/libint/$target_libint/$compiler/$compiler_version/nompi/include/libderiv
	cd -
  fi

# Preparation du fichier arch.popt 
  cd $compile_dir/$version/
  pwd
  echo " Preparation du fichier occigen.popt"

  if [ -e $compile_dir/$version/arch/$machine.popt ]; then
	echo " $machine.popt existe deja"
	set -x
	ls -la $machine.popt
 	rm -f $compile_dir/$version/arch/$machine.popt
  fi
  
  if [[ $version == *"2.5.1"* ]] ; then
	export LIBS_ADDS_ON=-L$compile_dir/$version/tools/hfx_tools/libint_tools/libint_cpp_wrapper.o
	else
	export LIBS_ADDS_ON="" 
  fi

# Utilisation de libint par standard_install.py 
  export LIBINT=/opt/software/$machine/libraries/libint/$target_libint/$compiler/$compiler_version
  #export DFLAGS="-D__INTEL -D__FFTSG -D__parallel  -D__MKL -D__FFTW3  -D__BLACS   -D__SCALAPACK   -D__LIBINT -D__LIBINT_MAX_AM=5 -D__HAS_NO_ISO_C_BINDING"
  export DFLAGS="-D__INTEL -D__INTEL_COMPILER  -D__parallel  -D__MKL -D__FFTW3  -D__BLACS   -D__SCALAPACK   -D__LIBINT -D__LIBINT_MAX_AM=5 -D__HAS_NO_ISO_C_BINDING"

  if [[ $libmpi == *"intelmpi"* ]] ; then
	export lmpi=$libmpi
	wrapper_compiler=mpiifort
	wrapper_LD=mpiifort
  else
	export lmpi=openmpi
	wrapper_compiler=mpif90
	wrapper_LD=mpif90
  fi

cat >arch_`date +"%m-%d-%y"`.popt <<EOF

CC       = icc
CXX      = icpc
CPP      = 
FC       = $wrapper_compiler -nofor-main
LD       = $wrapper_LD -nofor-main
AR       = ar -r

#CPPFLAGS = -C -traditional ${DFLAGS} -I$FFTW3_INC_DIR -I$LIBINT/include/libderiv  -I$LIBINT/include/libint  -I$LIBINT/include/libr12
CPPFLAGS = -traditional ${DFLAGS} -I$FFTW3_INC_DIR -I$LIBINT/include/libderiv  -I$LIBINT/include/libint  -I$LIBINT/include/libr12
#FCFLAGS  = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -I$MKL_INC_DIR/fftw -O2 -xCORE-AVX2 -heap-arrays 64 -funroll-loops -fpp -free 
FCFLAGS  = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -I$MKL_INC_DIR/fftw -O3 -xCORE-AVX2 -heap-arrays 64 -funroll-loops -fpp -free -fp-model precise -fp-model source -prec-div -prec-sqrt -fpe0
#FCFLAGS2 = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -I$MKL_INC_DIR/fftw -O1 -xCORE-AVX2 -heap-arrays 64 -fpp -free
FCFLAGS2 = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -I$MKL_INC_DIR/fftw -O1 -xCORE-AVX2 -heap-arrays 64 -fpp -free  
LDFLAGS  = ${FCFLAGS} -I$MKL_INC_DIR -nofor-main
#LIBS     =  $MKL_LIBS -lmkl_scalapack_lp64  -lmkl_blacs_${lmpi}_lp64    $FFTW3_LDFLAGS -L$LIBINT/lib -lderiv -lint -lr12  -lz -ldl -lstdc++   $LIBS_ADDS_ON
LIBS     = -nofor-main -L$MKL_LIB_DIR -lfftw3xf_intel $MKL_LIBS -lmkl_scalapack_lp64  -lmkl_blacs_${lmpi}_lp64     -L$LIBINT/lib -lderiv -lint -lr12  -lz -ldl -lstdc++   $LIBS_ADDS_ON
EOF
# With FCFLAGS (with O2) and FCFLAGS2 (with O1) lead to an instability in SCF calculation
#FCFLAGS  = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -O2 -xCORE-AVX2 -heap-arrays 64 -funroll-loops -fpp -free
#FCFLAGS2 = -debug -g ${DFLAGS} -I$MKL_INC_DIR  -O1 -xCORE-AVX2 -heap-arrays 64 -fpp -free
#FCFLAGS3 = -debug -g $(DFLAGS) -I$MKL_INC_DIR  -I$MKL_INC_DIR/fftw -O0 -traceback -fpp -free  

set +x
  more  arch_`date +"%m-%d-%y"`.popt 
  sleep 1
# Copie du nouveau fichier arch constitue dans le repertoire arch de la version
  ls arch_`date +"%m-%d-%y"`.popt
  rm -rf  $compile_dir/$version/arch/$machine.popt
  cp `pwd`/arch_`date +"%m-%d-%y"`.popt $compile_dir/$version/arch/$machine.popt	 
}

Install(){
  echo "************************************************************"
  echo "* Install"
  echo "************************************************************"
  cd $compile_dir/$version/makefiles
  pwd
  echo " Make clean"
  make -j24  ARCH=occigen VERSION=popt clean
  echo " Make "
  make -j4  ARCH=occigen VERSION=popt >& make.log
  cat make.log
set -x
# Liste des executables obtenus
  ls -lrt $compile_dir/$software_$version/exe/occigen

# ......................................
# Recopies dans le repertoire d'install
# ......................................
#for rep in arch data exe lib obj tests tools src
  for rep in arch data exe lib obj tests tools
  do
	if [ -d $install_dir/$rep ];then  
		rm -rf  $install_dir/$rep
    	fi
  	cp -rp $compile_dir/$software_$version/$rep $install_dir/.
  done
  echo " Contenu du repertoire d install exe"
  ls -la $install_dir/exe/occigen/.
  cp -r $install_dir/exe/occigen/ $install_dir/bin
  ls -lrt $install_dir/bin
  echo " Contenu du repertoire de tests"
  ls -la $install_dir/tests/.
  cd $my_sweet_home
}

#________________________________________________________________________
#________________________________________________________________________
Test(){
#...................................................
# Boucle sur quelques cas d'études de la methode QS
# indiques dans le fichier liste_test_cases
##...................................................
pwd
method=QS
  echo "************************************************************"
  echo "* Test sur $method "
  cp $my_sweet_home/liste_test_cases_input_files_QS $compile_dir
  echo "************************************************************"
#
for case_study in `cat liste_test_cases_input_files_$method` 
do
case="`echo $case_study | awk -F'/' {'print $1'}`"
#
#rm $case_`date +"%m-%d-%y"`.slurm 
#
tod="`date +"%m-%d-%y"`.slurm"
name="$case-$tod" 
echo $name
cat >$name <<EOF
#!/bin/bash
#SBATCH -J $scase
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=00:05:00
#SBATCH --constraint=HSW24
#SBATCH --output=$name.output

module load $compiler/$compiler_version
module load $libmpi/$mpi_release
module list

#ulimit -s unlimited
ulimit -s 256000
ulimit -c unlimited
ulimit -v unlimited
ulimit -a
srun --mpi=pmi2 -K1 --resv-ports -n \$SLURM_NTASKS $install_dir/exe/occigen/cp2k.popt $install_dir/tests/$method/$case_study
EOF
#
# Submit the job $case 
ls  -la $name 
sbatch  $name 
done
echo "Penser  A mettre A jour le script type donne par le module de cp2k" 
tail -20 *output
#---------------------------------------------------
}
#________________________________________________________________________
#________________________________________________________________________

RegTest(){
#.......................................................................................................
# CP2K comes with over 2500 test input files (located in cp2k/tests) which serve as both examples on how 
# to use the many features in CP2K and also as a method for developers to test modifications and 
# extensions to CP2K. In order to reduce the chance of bugs being introduced into the code, and ensure 
# that all parts of the code are working, we recommend that all users complete a test before using a 
# self-compiled binary for their projects 
#.......................................................................................................
	tod="`date +"%m-%d-%y"`.slurm"
	regTestName=regTest-popt_$tod

	if [[ $version == *"4.1"* ]] ; then
		CP2K_DIR=$compile_dir/cp2k_4.1/cp2k
	else
		CP2K_DIR=$compile_dir/cp2k-5.1
	fi

	cat > $regTestName <<EOF
#!/bin/bash
#SBATCH -J cp2k_reg_test_popt
#SBATCH --nodes=2
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=14
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --output regtest_cp2k_popt_%j.output
#SBATCH --constraint=BDW28
module purge
module load $compiler/$compiler_version
module load $libmpi/$mpi_release
module list

export CP2K_HOME=$CP2K_DIR
cd $CP2K_HOME
export PATH=$CP2K_HOME/exe/$machine:$PATH
export LD_LIBRARY_PATH=$CP2K_HOME/lib/$machine:$LD_LIBRARY_PATH
export LIBRARY_PATH=$CP2K_HOME/lib/$machine:$LIBRARY_PATH
export CP2K_DATA_DIR=/opt/software/occigen/applications/cp2k/4-DATA

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export cp2k_run_prefix="srun --mpi=pmi2 -K1  -m block:block -c $SLURM_CPUS_PER_TASK --resv-ports -n $SLURM_JOB_NUM_NODES"

time tools/regtesting/do_regtest -cp2kdir . -nosvn -nobuild -mpiranks $SLURM_JOB_NUM_NODES -ompthreads $OMP_NUM_THREADS  -maxtasks 16 -arch local -version popt
EOF

	sbatch $regTestName
}


Clean(){
  echo "************************************************************"
  echo "* clean : destruction du repertoire compile_dir"
  echo "************************************************************"
#  make clean
	rm -rf *xyz *wfn* *cube *restart *bak* *Hessian* *ener *RESTART *stress *slurm
	rm -rf $compile_dir

}
#------------------------------------------------------------------
#------------------------------------------------------------------
Deploy(){
  echo "************************************************************"
  echo "* Compilation de "$software" version "$version
  echo "*    compilateur "$compiler $compiler_version
  echo "*    libmpi      "$libmpi $mpi_release
  echo "************************************************************"
  islocal=`grep compilation config.cfg | awk -f'=' {'print $2'}`
#	 if [[ $islocal == *"local"* ]] ; then
#  		if [ -d $compile_dir ];then  
#    			rm -rf $compile_dir
#  		fi  
#		mkdir -p $compile_dir
#	fi

  if [ -d $install_dir ];then  
    rm -rf $install_dir
  fi  
  mkdir -p $install_dir

  if [ ! -d $compile_dir ];then
    mkdir -p $compile_dir
  fi
  
 Untar
 Configure_and_prepare
 Install
# Test
# Clean
}

Deploy
