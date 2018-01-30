#!/bin/bash

# Basically all this does is convert the lammps out file into a .xyz that is readable by VMD

sed "s/^1./C /g" -i position.xyz
sed "s/^2./O /g" -i position.xyz
sed "s/O 92/2592/g" -i position.xyz
