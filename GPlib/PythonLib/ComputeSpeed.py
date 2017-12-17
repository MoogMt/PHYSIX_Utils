# -*- coding: utf-8 -*-
"""
Computing VDOS

@author: CondensedOtters
"""

#============================
# Importing useful libraries
#==================================
import os
import platform
import numpy as np
import matplotlib.pyplot as plt
#==================================

#==================
# XYZ files
#========================================================
def countXYZstep( filepath_ , nb_atoms_ ):
    count = 0; 
    read = True;
    with open( filepath_, "r" ) as fp:
        while( read ): 
            for i in range(nb_atoms+2):
                line=fp.readline();
                if line == "":
                    read = False;
                    break;
                else: count += 1;
    return count/(nb_atoms+2);
def readXYZstep( file_pointer , nb_atoms , r_ ):
    for i in range(nb_atoms+2):
        line = file_pointer.readline();
        line_part = (line.rstrip("\n")).split()
        if line == "":
            return False
        if i >= 2:
            for j in range(3):
                r_[i-2,j] = line_part[j+1];
    return True
#========================================================

#========================
# Cell related functions
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
            dr[i,j] = minDir( r[i,j], r0[i,j] , cell[j] );
    return dr;
#========================================================

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

#=======================
# Physical parameters
#========================================================
# Number of dimension
ndim = 3;
# Step
step = 0;
# Timestep
timestep = 0.5;
# Sim print stride
sim_stride = 5;
# Timelaps
dt = timestep*sim_stride;
# Cell
a=9.0; b=9.0; c=9.0;
# Number of atoms
nb_atoms = 96;
# Nb of step in the simulations
nb_step = (int)(countXYZstep(filepath,nb_atoms));
# Reading Step parameters
start_step = 5000;
end_step = 10000000;
stride_comp = 1;
#--------------------------------------------------------
# Initiating variables of interest
#--------------------------------------------------------
# Names of atoms
name=[];
# Position at t
r  = np.zeros(( nb_atoms, 3 )); 
# Positions at t-dt
r0 = np.zeros(( nb_atoms, 3 )); 
# velocities x,y,z
v  = np.zeros(( nb_atoms, 3 ));
# storing velocities 
v_store = np.empty(( nb_atoms, 3, nb_step ));
#========================================================

#====================
# Reading TRAJEC.xyz
#========================================================
# Due to the repetition of lots of additions it might 
# actually be more efficient to first compute the number of 
# steps and initiate the v_store with that amount of memory 
# then fill it than dynamically filling it...
with open( filepath, "r" ) as fp:
    # Reading first step, initiates atomic positions
    if readXYZstep( fp, nb_atoms, r0 ) == False :
        print("Error Reading File!")
    # Reading all other steps  
    while( readXYZstep( fp, nb_atoms, r ) != False & step <= end_step ):
        # Compute speeds using finite elements
        if step > start_step:
            v = (minDist( r, r0, a, b, c )*1e-9)/(dt*1e-15);
            # Storing velocities in a vector
            v_store[:,:,step] = v
        # Remembers position for next step
        r0 = np.copy( r ); 
        # Incrementing steps
        print(step)
        step += 1;
#========================================================
print(np.max(v_store));
# Computing number of steps
nb_steps = v_store.size/v.size;

# DOING OPERATION ON THE VECTOR
from scipy.fftpack import fft, dct
from scipy import signal
    
vdos=np.zeros(nb_step*2-1);
for i in range(nb_atoms):
    for j in range(ndim):
        vdos = np.add(vdos,signal.correlate(v_store[i,j,:],v_store[i,j,:]))

c_atoms = np.arange(0,31,1)
vdos_C=np.zeros(nb_step*2-1);
for i in c_atoms:
    for j in range(ndim):
        vdos_C = np.add(vdos_C,signal.correlate(v_store[i,j,:],v_store[i,j,:]))
        
o_atoms = np.arange(32,95,1)
vdos_O=np.zeros(nb_step*2-1);
for i in o_atoms:
    for j in range(ndim):
        vdos_O = np.add(vdos_O,signal.correlate(v_store[i,j,:],v_store[i,j,:]))
        
x = np.arange(0, nb_step*2-1, 1);
plt.xlim([0,45000])
plt.plot(x,dct(vdos))
plt.figure();
plt.plot(x,dct(vdos_C))
plt.figure();
plt.plot(x,dct(vdos_O))


