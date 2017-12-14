# -*- coding: utf-8 -*-
"""
COMPUTING SPEED

@author: CondensedOtters


"""

# Importing useful libraries
import os
import numpy as np
import matplotlib.pyplot as plt

# Folder parameters
#====================================================================
# Folder containing all files
folder="/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2000K/"
# Local File
filepath=folder+"TRAJEC_wrapped.xyz"
#====================================================================

# Physical parameters
#========================================================
# Step
step = 0;
# Timestep
timestep = 0.5;
# Timelaps
timelaps = timestep*5;
# Cell
a=9.0; b=9.0; c=9.0;
# Number of atoms
nb_atoms = 96;
#========================================================

#==================
# Reading XYZ step
#========================================================
def readXYZstep( file_pointer , nb_atoms , x_ , y_ , z_):
    for i in range(nb_atoms+2):
        line = file_pointer.readline()
        line_part = (line.rstrip("\n")).split()
        if line == "":
            return False
        if i >= 2:
            x_[i-2] = line_part[1];
            y_[i-2] = line_part[2];
            z_[i-2] = line_part[3];
    return x_, y_ ,z_ 
#========================================================

#============
# Atom Names
#========================================================
name=[];
# Position x,y,z
x=np.zeros(nb_atoms); x0=np.zeros(nb_atoms); 
y=np.zeros(nb_atoms); y0=np.zeros(nb_atoms);
z=np.zeros(nb_atoms); z0=np.zeros(nb_atoms);
# velocities x,y,z
vx=np.array([]);
vy=np.array([]);
vz=np.array([]);
#========================================================

#====================
# Reading TRAJEC.xyz
#========================================================
count=0;
with open(filepath,"r") as fp:
    # Reading first step
    if readXYZstep(fp,nb_atoms,x0,y0,z0) == False :
        print("Error Reading File!")
    # Reading all other steps  
    while( readXYZstep(fp,nb_atoms,x,y,z) != False ):
        dx=x0[0]-x[0];
        dy=y0[0]-y[0];
        dz=z0[0]-z[0];
        x0=np.copy(x); 
        y0=np.copy(y);
        z0=np.copy(z);
        print(step)
        step+=1;
#========================================================