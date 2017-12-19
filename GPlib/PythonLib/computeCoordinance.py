# -*- coding: utf-8 -*-
"""
Computing VDOS

@author: CondensedOtters
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
#--------------------------------------------------------
#--------------------------------------------------------
def minDist( r, r0, cell):
    if type(r) == float :
        dr = 0;
        for j in range( len( cell ) ):
            dr[j] = minDir( r[j], r0[j] , cell[j] );
    else:
        dr = np.zeros(( r[:,0].size, 3 ));
        for i in range(r[:,0].size):
            for j in range( len( cell )-3 ):
                dr[i,j] = minDir( r[i,j], r0[i,j] , cell[j] );
    return dr;
#========================================================

#=======================
# Physical parameters
#========================================================
ndim = 3; 
femto = 1e-15;
pico  = 1e-12;
angstrom = 1e-10;
nano  = 1e-9;
micro = 1e-6;
milli = 1e-3;
kilo  = 1e3;
mega  = 1e6;
giga  = 1e9;
tera  = 1e12
peta  = 1e15;
#====================================================================

#===========================
# Sigmoid PLUMED function
#===================================================================
def sigm(x_,x0_,n_,m_):
    return (1-(x_/x0_)**n_)/(1-(x_/x0_)**m_);
def sigMatrix(matrix_,x0_,n_,m_):
    if type( matrix_ ) == int | type(matrix_) == float :
        return sigm(matrix_,x0_,n_,m_);
    else:
        sig=np.vectorize(sigm);
        return sig(matrix_,x0_,n_,m_)
#====================================================================



#===========================================================================
def computeCoordiance(r_, cell_ ):
    # Sigmoid parameters
    x0_ =  1.8; n_ = 8; m_= 16;
    return sigMatrix(computeContactMatrix(r_),x0_,n_,m_);
#----------------------------------------------------------------------            
def computeCoord( filepath_ , nb_atoms_ , start_step_, end_step_, cell_):
    #================
    # Init variables
    #========================================================
    step_ = 0; # Number of step of the simulation
    nb_step_ = (int)( countXYZstep( filepath_, nb_atoms_ ) ); 
    r  = np.zeros(( nb_atoms_, 3 ));
    # Coordinances
    coord_avg_  = np.zeros( nb_step_ );
    coord_list_ = np.zeros(( nb_atoms_, nb_step_ )); 
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    # Opening file
    with open( filepath_, "r" ) as fp_:
        # Reading file
        while( readXYZstep( fp_, nb_atoms_, r ) != False & step_ <= end_step_ ):
            # Compute speeds using finite elements
            if step_ >= start_step_:
                #DO STUFF HERE
                coord_list_[:,step_-start_step_] = computeCoordiance(r);
                coord_avg_[step_-start_step_] = np.mean(coord_list_[:,step_-start_step_]);
            # Incrementing steps
            print(step_)
            step_ += 1;
    #========================================================
    
    return coord_avg_, coord_list_;
#====================================================================





