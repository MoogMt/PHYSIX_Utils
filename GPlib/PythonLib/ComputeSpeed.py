# -*- coding: utf-8 -*-
"""
COMPUTING SPEED

@author: CondensedOtters


"""

# Importing useful libraries
import os
import platform
import numpy as np
import matplotlib.pyplot as plt

# Determining P,T to target and if sprint
#====================================================================
T = "2000";
Length_cell = "9.0";
SPRINT = False;
#====================================================================1

# Folder parameters
#====================================================================
# Rough attempt at automating target file detection regardless of OS
folder=""
if platform.system() == "Windows":
    folder = "E:\Data\CO2\AIMD\Liquid\PBE-MT\9.0\\2000K\\"
elif platform.system() == "Linux":
    username = os.getusername()
    path_in_usb = "/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2000K"
    folder = "/media/" + username + path_in_usb;
# Targetting file
file = "TRAJEC_wrapped.xyz"
# Joining folder and file
filepath = os.path.join(folder,file)
#====================================================================

# Physical parameters
#========================================================
# Step
step = 0;
# Timestep
timestep = 0.5;
# Timelaps
dt = timestep*5;
# Cell
a=9.0; b=9.0; c=9.0;
# Number of atoms
nb_atoms = 96;
# Names of atoms
name=[];
# Position at t
r  = np.zeros(( nb_atoms, 3 )); 
# Positions at t-dt
r0 = np.zeros(( nb_atoms, 3 )); 
# velocities x,y,z
v  = np.zeros(( nb_atoms, 3 ));
# storing velocities 
v_store = np.empty(( nb_atoms, 3 ));
#========================================================
#========================================================

#==================
# Reading XYZ step
#========================================================
def readXYZstep( file_pointer , nb_atoms , r_ ):
    for i in range(nb_atoms+2):
        line = file_pointer.readline()
        line_part = (line.rstrip("\n")).split()
        if line == "":
            return False
        if i >= 2:
            for j in range(3):
                r_[i-2,j] = line_part[j+1];
    return True
#========================================================

#========================================================
def minDir( x , x0, a ):
    dx=x-x0;
    if dx > a*0.5: 
        return dx-a
    elif dx<-a*0.5: 
        return dx+a;
    else: 
        return dx;
#--------------------------------------------------------
def minDist( r, r0, a, b, c ):
    dr = np.zeros(( r[:,0].size, 3 ));
    cell=[ a, b, c ];
    for i in range(r[:,0].size):
        for j in range( len( cell ) ):
            dr[i,j] = minDir( r[i,j], r0[i,j] , cell[j] )
    return dr;
#========================================================

#====================
# Reading TRAJEC.xyz
#========================================================
with open( filepath, "r" ) as fp:
    # Reading first step
    if readXYZstep( fp, nb_atoms, r0 ) == False :
        print("Error Reading File!")
    # Reading all other steps  
    while( readXYZstep( fp, nb_atoms, r ) != False ):
        # Compute speeds using finite elements
        v = (minDist( r, r0, a, b, c ))/dt;
        if step < 10:
            np.append( v_store, v);
        # Remembers last positions
        r0=np.copy( r ); 
        # Incrementing steps
        print(step)
        step += 1;
print(v_store.size);
#========================================================
        
