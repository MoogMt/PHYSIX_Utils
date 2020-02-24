#!/bin/bash

for volume in 8.82 9.0 9.05 9.1 9.2 9.3 9.35 9.375 9.4 9.5 9.8
do
    cd $volume"/3000K/"
    cp ../../piv_clustering.x  .
    mpirun -np 4 ./piv_clustering.x -filexyz traj_stride.xyz -bsize $volume $volume $volume -method 1 -coord1_range 1.2 2.5 -algorithm 20
    cd ../..
done
