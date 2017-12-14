# -*- coding: utf-8 -*-
"""
COMPUTING SPEED

@author: CondensedOtters


"""

# Importing useful libraries
import numpy as np

# Folder containing all files
folder="/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/"

# Local File
filepath=folder+"TRAJEC.xyz"

#========================================================
# Atom Names
name=[];
# Position x,y,z
x=np.array([]), y=np.array([]), z=np.array([]);
# velocities x,y,z
vx=np.array([]), vy=np.array([]), vz=np.array([]);
#========================================================

# Reading TRAJEC.xyz
with open(filepath) as fp:
    # Getting the line
    line = fp.readline();
    # Removing the return
    line = line.rsplit('\n')
    # Parsing names, x, y, z
    
    np.append(x,line[1]);
    np.append(y,line[2]);
    np.append(z,line[3]);
