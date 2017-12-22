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
def minDir( x_ , x0_, a_ ):
    dx_ = x_ - x0_;
    if dx_ > a_*0.5: 
        return dx_-a_
    elif dx_< - a_*0.5: 
        return dx_ + a_;
    else: 
        return dx_;
#========================================================
        
#========================================================
def addVector( r_, r0_ ):
    dr_ = np.zeros(( r_[:,0].size, 3 ));
    for i in range(r_[:,0].size):
        for j in range( 3 ):
            dr_[i,j] = r_[i,j] + r0_[i,j];
    return dr_;
#========================================================

#=======================
# Physical parameters
#========================================================
femto = 1e-15;
Tera = 1e12
Thz2cm = 100/2.99792;
angstrom = 1e-10;
#====================================================================

#===================================================================
def computeVDOS(filepath_, nb_atoms_, start_step_, end_step_, dt_ ):
    # Step
    step_ = 0;
    # Nb of step in the simulations
    nb_step_ = (int)( countXYZstep( filepath_, nb_atoms_ ) );
    # Position at t
    r  = np.zeros(( nb_atoms_, 3 )); 
    # Positions at t-dt
    r0 = np.zeros(( nb_atoms_, 3 )); 
    # velocities x,y,z
    v  = np.zeros(( nb_atoms_, 3 ));
    # storing velocities 
    v_store = np.zeros(( nb_atoms_, 3, (nb_step_-start_step_)-1 ));
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    with open( filepath_, "r" ) as fp:
        # Reading first step, initiates atomic positions
        if readXYZstep( fp, nb_atoms_, r0 ) == False :
            print("Error Reading File!")
            # Reading all other steps  
        while( readXYZstep( fp, nb_atoms_, r ) != False & step_ <= end_step_ ):
            # Compute speeds using finite elements
            if step_ >= start_step_:
                v = addVector( r, -r0 )/dt_*angstrom/femto;
                # Storing velocities in a vector
                v_store[:,:,step_-start_step_] = v
            # Remembers position for next step
            r0 = np.copy( r ); 
            # Incrementing steps
            print(step_)
            step_ += 1;
    #========================================================
    
    #================
    # Computing VDOS
    #========================================================
    vdos_=np.zeros((nb_step_-start_step_-1));
    for i in range(nb_atoms_):
        v_ = v_store[i,0,0]*v_store[i,0,:]+v_store[i,1,0]*v_store[i,1,:]+v_store[i,2,0]*v_store[i,2,:];
        vdos_ = np.add(vdos_,signal.correlate(v_,v_,mode="same",method="fft"))
    # Normalization + FFT (DCT type 2)
    vdos_ = abs(dct(vdos_/(nb_atoms_),type=2,norm="ortho"));
    return  nb_step_-start_step_, vdos_
#========================================================
    
#========================
# Simulation parameters
#========================================================
nb_atoms = 96;          # Number of atoms
start_step = 5000;      # Start step
end_step = 1000000000;  # End step
timestep = 0.483776856; # CPMD timestep
sim_stride = 5;         # Stride of printing position in simulation
dt=5*0.483776856;       # efficient dt
a_882 = 8.82; a_900 = 9.00; a_910 = 9.10; a_980 = 9.80 # Cell parameters
#========================================================

#=======================
# Computing VDOS vs P
#========================================================
nbstep_882, vdos_882 = computeVDOS("/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2000K/TRAJEC.xyz", nb_atoms, start_step, end_step, dt );
nbstep_900, vdos_900 = computeVDOS("/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2000K/TRAJEC.xyz",  nb_atoms, start_step, end_step, dt );
nbstep_910, vdos_910 = computeVDOS("/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.1/2000K/TRAJEC.xyz",  nb_atoms, start_step, end_step, dt );
nbstep_980, vdos_980 = computeVDOS("/media/moogmt/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.8/2000K/TRAJEC.xyz",  nb_atoms, start_step, end_step, dt );
#========================================================

#====================
# Plotting VDOS vs P
#========================================================
x882 = np.arange(0, vdos_882.size , 1);
x900 = np.arange(0, vdos_900.size , 1);
x910 = np.arange(0, vdos_910.size , 1);
x980 = np.arange(0, vdos_980.size , 1);
plt.plot((x980/(nbstep_980*dt*femto))/2/Tera*Thz2cm,vdos_980+140000,'c.');
plt.plot((x910/(nbstep_910*dt*femto))/10/Tera*Thz2cm,vdos_910*5+100000,'g.');
plt.plot((x900/(nbstep_900*dt*femto))/2/Tera*Thz2cm,vdos_900+50000,'b.');
plt.plot((x882/(nbstep_882*dt*femto))/2/Tera*Thz2cm,vdos_882,'r.');
plt.legend(["40GPa - 2000K","50GPa - 2000K (9.1)", "50GPa - 2000K (9.0)","60GPa - 2000K"])
plt.ylabel("VDOS (arb. u.)");
plt.xlabel("Wavenumber (cm-1)");
plt.xlim([0,3000])
plt.show();
#========================================================

