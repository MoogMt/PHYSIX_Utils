#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:41:31 2017

@author: moogmt
"""

#============================
# Importing useful libraries
#==================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
from scipy import signal
#==================================

#==================
# XYZ files
#========================================================
def countXYZstep( filepath_ , nb_atoms_ ):
    count = 0; 
    read = True;
    with open( filepath_, "r" ) as fp:
        while( read ): 
            for i in range(nb_atoms_+2):
                line=fp.readline();
                if line == "":
                    read = False;
                    break;
                else: count += 1;
    return count/(nb_atoms_+2);
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
#========================================================

#=======================
# Physical parameters
#========================================================
# Number of dimension
ndim = 3;
femto = 1e-15;
Tera = 1e12
Thz2cm = 100/2.99792;
angstrom = 1e-10;
#====================================================================

def computeVelocities( filepath_, nb_atoms_, start_step_, end_step_, dt_ ):
    # Step
    step_ = 0;
    # Nb of step in the simulations
    nb_step_ = (int)( countXYZstep( filepath_, nb_atoms_ ) );
    # Position at t
    r_  = np.zeros(( nb_atoms_, 3 )); 
    # Positions at t-dt
    r0_ = np.zeros(( nb_atoms_, 3 )); 
    # velocities x,y,z
    v_  = np.zeros(( nb_atoms_, 3 ));
    # storing velocities 
    v_store = np.zeros(( nb_atoms_, 3, (nb_step_-start_step_)-1 ));
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    with open( filepath_, "r" ) as fp:
        # Reading first step, initiates atomic positions
        if readXYZstep( fp, nb_atoms_, r0_ ) == False :
            print("Error Reading File!")
            # Reading all other steps  
        while( readXYZstep( fp, nb_atoms_, r_ ) != False & step_ <= end_step_ ):
            # Compute speeds using finite elements
            if step_ >= start_step_:
                v_ = np.add( r_, -r0_)*angstrom/(dt_*femto);
                # Storing velocities in a vector
                v_store[:,:,step_-start_step_] = v_
            # Remembers position for next step
            r0_ = np.copy( r_ ); 
            # Incrementing steps
            print(step_)
            step_ += 1;
    #========================================================
    
    return v_store

def writeVelocities( filepath_ , v_store):
    truc = open( filepath_,"w")
    for i in range(v_store[0,0,:].size):
        for j in range(v_store[:,0,0].size):
            for k in range(v_store[0,:,0].size):
                truc.write(str(v_store[j,k,i]));
                truc.write(" ");
            truc.write("\n")
        truc.write(v_store[0,0,:]);
        truc.write("\n")
    truc.close() 
    return

nb_atoms = 96;
start_step = 2000;
end_step   = 10000000000000
dt = 0.48*5;

writeVelocities( "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2000K/velocities.dat", computeVelocities(  "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2000K/TRAJEC.xyz", nb_atoms, start_step, end_step, dt ))
