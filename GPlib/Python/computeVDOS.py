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
 