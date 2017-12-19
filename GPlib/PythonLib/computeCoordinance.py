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
    if np.ndim(r) == 1 :
        dr = 0;
        for j in range( len( cell ) ):
            dr += minDir( r[j], r0[j] , cell[j] )**2;
        dr = np.sqrt(dr);
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
    return (1-(x_/x0_+0.0000001)**n_)/(1-(x_/x0_+0.0000001)**m_);
def sigMatrix(matrix_,x0_,n_,m_):
    if type( matrix_ ) == int or  type(matrix_) == float :
        return sigm(matrix_,x0_,n_,m_);
    else:
        sig=np.vectorize(sigm);
        return sig(matrix_,x0_,n_,m_)
#====================================================================

# Contact Matrix
#====================================================================
def computeContactMatrix( r_, cell_ ):
    x0_ =  1.8; n_ = 6; m_= 36;
    matrix_ = np.zeros((len(r_),len(r_)));
    for i in range( len(r_) ):
        for j in np.arange( i, len(r_), 1 ):
            value_ = minDist( r_[i,:], r_[j,:], cell_ );
            matrix_[i,j] = sigm( value_, x0_, n_, m_ ); 
            matrix_[j,i] = sigm( value_, x0_, n_, m_ );
    return matrix_;
#====================================================================

#==============================================================================
def computeCoordFromCM( matrix_ ):           
    coord_vec = np.zeros(( len( matrix_ ) ));
    for i in range( coord_vec.size ):
        value = 0;
        for j in range( coord_vec.size ):
            value += matrix_[i,j];
        coord_vec[i] = value;
    return 
#------------------------------------------------------------------------------
def computeCoord( filepath_ , nb_atoms_ , start_step_, end_step_, cell_ ):
    #================
    # Init variables
    #========================================================
    step_ = 0; # Number of step of the simulation
    nb_step_ = (int)( countXYZstep( filepath_, nb_atoms_ ) ); 
    r_  = np.zeros(( nb_atoms_, 3 ));
    # Coordinances
    coord_avg_  = np.zeros( nb_step_ );
    coord_list_ = np.zeros(( nb_atoms_, nb_step_-start_step )); 
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    # Opening file
    with open( filepath_, "r" ) as fp_:
        # Reading file
        while( readXYZstep( fp_, nb_atoms_, r_ ) != False & step_ <= end_step_ ):
            # Compute speeds using finite elements
            if step_ >= start_step_:
                #DO STUFF HERE
                contact_matrix_ = computeContactMatrix(r_,cell_);
                coord_list_[:,step_-start_step_] = computeCoordFromCM( contact_matrix_ );
                coord_avg_[step_-start_step_] = np.mean(coord_list_[:,step_-start_step_]);
            # Incrementing steps
            print(step_)
            step_ += 1;
    #========================================================
    
    return coord_avg_, coord_list_;
#====================================================================

nb_atoms = 96;
start_step = 5000;
end_step = 10000000;
cell = np.array([8.82,8.82,8.82])

coord882 = computeCoord( "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2000K/TRAJEC_wrapped.xyz" , nb_atoms , start_step, end_step, cell)


