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
def minDist( r, r0, a, b, c ):
    dr = np.zeros(( r[:,0].size, 3 ));
    cell=[ a, b, c ];
    for i in range(r[:,0].size):
        for j in range( len( cell ) ):
            dr[i,j] = minDir( r[i,j], r0[i,j] , cell[j] );
    return dr;
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

#===================================================================
def computeVDOS(filepath, nb_atoms_, a , b, c, ndim, start_step, end_step, timestep, sim_stride ):
    # Step
    step = 0;
    # Nb of step in the simulations
    nb_step = (int)(countXYZstep(filepath,nb_atoms_));
    #--------------------------------------------------------
    # Initiating variables of interest
    #--------------------------------------------------------
    # Position at t
    r  = np.zeros(( nb_atoms_, 3 )); 
    # Positions at t-dt
    r0 = np.zeros(( nb_atoms_, 3 )); 
    # velocities x,y,z
    v  = np.zeros(( nb_atoms_, 3 ));
    # storing velocities 
    v_store = np.zeros(( nb_atoms_, 3, (nb_step-start_step)-1 ));
    # 
    dt = timestep*sim_stride;
    #========================================================
    
    #====================
    # Reading TRAJEC.xyz
    #========================================================
    with open( filepath, "r" ) as fp:
        # Reading first step, initiates atomic positions
        if readXYZstep( fp, nb_atoms_, r0 ) == False :
            print("Error Reading File!")
            # Reading all other steps  
        while( readXYZstep( fp, nb_atoms_, r ) != False & step <= end_step ):
            # Compute speeds using finite elements
            if step >= start_step:
                v = minDist( r, r0, a, b, c )/dt*angstrom/femto;
                # Storing velocities in a vector
                v_store[:,:,step-start_step] = v
            # Remembers position for next step
            r0 = np.copy( r ); 
            # Incrementing steps
            print(step)
            step += 1;
    #========================================================
    
    #================
    # Computing VDOS
    #========================================================
    vdos=np.zeros((nb_step-start_step-1));
    for i in range(nb_atoms_):
        for j in range(ndim):
            vdos = np.add(vdos,signal.correlate(v_store[i,j,:],v_store[i,j,:],mode="same",method="fft"))
            vdos = abs(dct(vdos/(ndim*nb_atoms_),type=2,norm="ortho"));
    vdos2=np.empty([])
    for k in range(vdos.size):
        if k%2 == 0:
            vdos2 = np.append(vdos2,vdos[k])
    vdos3=np.empty([])
    for i in range(vdos.size):
        if i%2 != 0:
            vdos3 = np.append(vdos2,vdos[i])
    return  nb_step-start_step, vdos2/(nb_step-start_step), vdos3/(nb_step-start_step)
#========================================================

""" Overly clever stuff
# Folder parameters
#====================================================================
# Rough attempt at automating target file detection regardless of OS
folder=""
if platform.system() == "Windows":
    folder = "E:\Data\CO2\AIMD\Liquid\PBE-MT\8.82\2500K\\"
elif platform.system() == "Linux":
    import getpass
    username = getpass.getuser()
    path_in_usb = "/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2500K"
    folder = "/media/" + username + path_in_usb;
# Targetting file
file = "TRAJEC_wrapped.xyz"
# Joining folder and file
filepath = os.path.join(folder,file)
#====================================================================
"""

dt=5*0.483776856

# Plotting vdos
#========================================================
nbstep_882, vdos_882, vdos= computeVDOS("/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/8.82/2500K/TRAJEC.xyz", 96, 8.82, 8.82, 8.82, ndim, 2000, 1000000000, 0.483776856, 5);
nbstep_900, vdos_900, vdos2 = computeVDOS("/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.0/2500K/TRAJEC.xyz",  96, 9.00, 9.00, 9.00, ndim, 2000, 1000000000, 0.483776856, 5);
#nbstep_910, vdos_910 = computeVDOS("/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.1/2000K/TRAJEC.xyz",  96, 9.10, 9.10, 9.10, ndim, 2000, 1000000000, 0.483776856, 5);
nbstep_980, vdos_980, vdos3 = computeVDOS("/media/moog/KINGSTON/Data/CO2/AIMD/Liquid/PBE-MT/9.8/2500K/TRAJEC.xyz",  96, 9.80, 9.80, 9.80, ndim, 2000, 1000000000, 0.483776856, 5);
x882 = np.arange(0, vdos_882.size , 1);
x900 = np.arange(0, vdos_900.size , 1);
#x910 = np.arange(0, vdos_910.size , 1);
x980 = np.arange(0, vdos_980.size , 1);
plt.ylabel("VDOS (arb. u.)");
plt.xlabel("Wavenumber (cm-1)");
plt.xlim([0,3000]);
plt.plot((x882/(nbstep_882*dt*femto))/Tera*Thz2cm,vdos_882+2000,'r-');
plt.plot((x900/(nbstep_900*dt*femto))/Tera*Thz2cm,vdos_900+1000,'b-');
#plt.plot((x910/(nbstep_910*dt*femto))/Tera*Thz2cm,vdos_910,'g.');
plt.plot((x980/(nbstep_980*dt*femto))/Tera*Thz2cm,vdos_980,'c-');
plt.legend(["60GPa - 2500K","50GPa - 2500K", "40GPa - 2500K"])
plt.show();

x882 = np.arange(0, vdos.size , 1);
x900 = np.arange(0, vdos2.size , 1);
#x910 = np.arange(0, vdos_910.size , 1);
x980 = np.arange(0, vdos3.size , 1);
plt.ylabel("VDOS (arb. u.)");
plt.xlabel("Wavenumber (cm-1)");
plt.xlim([0,3000]);
plt.plot((x882/(nbstep_882*dt*femto))/Tera*Thz2cm,vdos+2000,'r-');
plt.plot((x900/(nbstep_900*dt*femto))/Tera*Thz2cm,vdos2+1000,'b-');
#plt.plot((x910/(nbstep_910*dt*femto))/Tera*Thz2cm,vdos_910,'g.');
plt.plot((x980/(nbstep_980*dt*femto))/Tera*Thz2cm,vdos3,'c-');
plt.legend(["60GPa - 2500K","50GPa - 2500K", "40GPa - 2500K"])
plt.show();

#========================================================