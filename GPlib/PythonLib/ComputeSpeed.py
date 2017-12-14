# -*- coding: utf-8 -*-
"""
COMPUTING SPEED

@author: CondensedOtters


"""

# Importing useful libraries
import numpy as np

# Folder parameters
#====================================================================
# Folder containing all files
folder="/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2000K/"
# Local File
filepath=folder+"TRAJEC.xyz"
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
#========================================================

#========================================================
# Atom Names
name=[];
# Position x,y,z
x=np.array([]); x0=np.array([]);
y=np.array([]); y0=np.array([]);
z=np.array([]); z0=np.array([]);
# velocities x,y,z
vx=np.array([]);
vy=np.array([]);
vz=np.array([]);
#========================================================

#====================
# Reading TRAJEC.xyz
#=================================
count=0;
with open(filepath,"r") as fp:
    # Reading a file line by line
    for line in fp: 
        if count >= 2:    
            # Parsing names, x, y, z
            line_split = line.split();
            name.append(line_split[0]);
            np.append(x,line_split[1]);
            np.append(x,line_split[2]);
            np.append(x,line_split[3]);
        if count == 96:
            print(step)
            step += 1;
            np.append(vx,(x-x0)/timelaps);
            np.append(vy,(y-y0)/timelaps);
            np.append(vz,(z-z0)/timelaps);
            x=x0; y=y0; z=z0;
            np.empty(x); np.empty(y); np.empty(z);
            x=x0; y=y0; z=z0;
            count=-2;
        count += 1;
#=================================