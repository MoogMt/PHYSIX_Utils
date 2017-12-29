# -*- coding: utf-8 -*-
"""
Computing VDOS

@author: CondensedOtters
"""

#============================
# Importing useful libraries
#==================================
import numpy as np
import numba as nb
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
def minDir( x_ , x0_, a_ ):
    dx_=x_-x0_;
    if dx_ > a_*0.5: 
        return dx_-a_
    elif dx_<-a_*0.5: 
        return dx_+a_;
    else: 
        return dx_;
#--------------------------------------------------------
#--------------------------------------------------------
def minDist( r_, r0_, cell_):
    if np.ndim(r_) == 1 :
        dr_ = 0;
        for j in range( len( cell_ ) - 3 ):
            dr_ += minDir( r_[j], r0_[j] , cell_[j] )**2;
        dr_ = np.sqrt(dr_);
    else:
        dr_ = np.zeros(( r_[:,0].size, 3 ));
        for i in range(r_[:,0].size):
            for j in range( len( cell_ )-3 ):
                dr_[i,j] = minDir( r_[i,j], r0_[i,j] , cell_[j] );
    return dr_;
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
        for j in np.arange( i+1, len(r_), 1 ):
            value_ = minDist( r_[i,:], r_[j,:], cell_ );
            matrix_[i,j] = sigm( value_, x0_, n_, m_ ); 
            matrix_[j,i] = sigm( value_, x0_, n_, m_ );
    return matrix_;
#====================================================================

#==============================================================================
def computeCoordFromCM( matrix_ ):           
    coord_vec_ = np.zeros(( len( matrix_ ) ));
    for i in range( coord_vec_.size ):
        value = 0;
        for j in range( coord_vec_.size ):
            value += matrix_[i,j];
        coord_vec_[i] = value;
    return coord_vec_
#------------------------------------------------------------------------------
def computeCoord( filepath_ , nb_atoms_ , start_step_, end_step_, cell_ ):
    #================
    # Init variables
    #========================================================
    step_ = 0; # Number of step of the simulation
    nb_step_ = (int)( countXYZstep( filepath_, nb_atoms_ ) ); 
    r_  = np.zeros(( nb_atoms_, 3 ));
    # Coordinances
    coord_list_ = np.zeros(( nb_atoms_, nb_step_- start_step_ )); 
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    # Opening file
    with open( filepath_, "r" ) as fp_:
        # Reading file
        while( readXYZstep( fp_, nb_atoms_, r_ ) != False & step_ <= end_step_ ):
            # Compute speeds using finite elements
            if step_ >= start_step_ :
                contact_matrix_ = computeContactMatrix(r_,cell_);
                coord_list_[:,step_-start_step_] = computeCoordFromCM( contact_matrix_ ); 
            # Incrementing steps
            print(step_)
            step_ += 1;
    #========================================================
    
    return coord_list_;
#====================================================================

#-----------------
# Local paramters
#----------------------------------------------
nb_atoms = 96;
start_step = 5000;
end_step = 100000000;
sim_stride = 5;
sim_dt = 0.00048;
dt = sim_stride*sim_dt;
#----------------------------------------------
# Cells
cell882 = np.array([8.82,8.82,8.82,90,90,90]);
cell900 = np.array([9,9,9,90,90,90]);
cell980 = np.array([9.8,9.8,9.8,90,90,90]);
#----------------------------------------------

#------------------------
# Computing coordinance
#----------------------------------------------
print("Now Dealing with 8.82 - 2000K");
coord_list_882 = computeCoord( "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2000K/TRAJEC_wrapped.xyz" , nb_atoms , start_step, end_step, cell882)
print("Now Dealing with 9.00 - 2000K");
coord_list_900 = computeCoord( "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2000K/TRAJEC_wrapped.xyz" , nb_atoms , start_step, end_step, cell900)
print("Now Dealing with 9.80 - 2000K");
coord_list_980 = computeCoord( "/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.8/2000K/TRAJEC_wrapped.xyz" , nb_atoms , start_step, end_step, cell980)
#----------------------------------------------

coord_total_882 = coord_list_882.reshape(coord_list_882[:,:].size)
coord_total_900 = coord_list_900.reshape(coord_list_900[:,:].size)
coord_total_980 = coord_list_980.reshape(coord_list_980[:,:].size)

hist882 = np.histogram(coord_total_882,bins=200,normed=True)[1]
hist900 = np.histogram(coord_total_900,bins=200,normed=True)[1]
hist980 = np.histogram(coord_total_980,bins=200,normed=True)[1]

plt.hist(hist882)
plt.hist(hist900)
plt.hist(hist980)