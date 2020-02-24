#!/bin/bash

#SBATCH -J relax
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=02:00:00
#SBATCH --output relax.out

module purge
module load intel/15.0.3.187
module load intelmpi/5.0.1.035

EXE=patate
SEED=loutre

srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS $EXE $SEED

rm $SEED.castep_bin $SEED.bands $SEED.check $SEED.geom $SEED.cst_esp
