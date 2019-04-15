#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 22:04:24 2019

@author: moogmt
"""

import numpy as np
import matplotlib.pyplot as plt
   

#=======================
# Physical parameters
#========================================================
femto = 1e-15;
Tera = 1e12
Thz2cm = 100/2.99792;
angstrom = 1e-10;
#====================================================================

#========================
# Simulation parameters
#========================================================
nb_atoms = 96;          # Number of atoms
start_step = 0;      # Start step
end_step = 20000;  # End step
timestep = 2*0.483776856; # CPMD timestep
sim_stride = 5;         # Stride of printing position in simulation
dt=sim_stride*timestep;       # efficient dt
a=9.8
#========================================================

#=======================
# Computing VDOS vs P
#========================================================
nbstep, vdos = computeVDOS("/home/moogmt/CO2/CO2_AIMD/9.8/3000K/TRAJEC_wrapped.xyz", nb_atoms, start_step, end_step, dt );

#====================
# Plotting VDOS vs P
#========================================================
cm98 = np.arange(0, vdos.size , 1);
plt.plot((cm98/(nbstep*dt*femto))/2/Tera*Thz2cm,vdos,'c.');
plt.legend(["20GPa - 2000K"])
plt.ylabel("VDOS (arb. u.)");
plt.xlabel("Wavenumber (cm-1)");
plt.xlim([0,3000])
plt.show();
#========================================================

