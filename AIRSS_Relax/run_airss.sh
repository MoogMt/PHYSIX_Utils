#!/bin/bash

#SBATCH -J AIRSS
#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --output AIRSS.out

module purge
module load intel/15.0.3.187
module load intelmpi/5.0.1.035

PRES=30
SEED=4CO2
cmd="airss.pl -mpinp 48 -press "$PRES" -seed "$SEED" < /dev/null >&/dev/null &"
$cmd

