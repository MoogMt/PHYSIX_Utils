#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:51:00 2019

@author: moogmt
"""


import matplotlib.pyplot as plt #1
import numpy as np
from scipy.optimize import curve_fit

end_P_1=47
end_P_2=70



P_fic1=np.linspace(25,end_P_1,100)
P_fic2=np.linspace(end_P_1,end_P_2,100)

folder=str("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/")

P_2000,D_2000 = np.loadtxt(str(folder+"MSD_PV-2000K.dat"), unpack=True )
P_2500,D_2500 = np.loadtxt(str(folder+"MSD_PV-2500K.dat"), unpack=True )
P_3000,D_3000 = np.loadtxt(str(folder+"MSD_PV-3000K.dat"), unpack=True )

plt.plot(P_2000,D_2000,"r-o")
plt.plot(P_2500,D_2500,"b-o")
plt.plot(P_3000,D_3000,"g-o")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()

def expdec(x,A,B,C):
    return B*np.exp(-A*(x-C))

#==============================================================================

end_P_1_3K=44

P_3000_1=P_3000[P_3000<end_P_1]
D_3000_1=D_3000[P_3000<end_P_1]

P_fic1_3K=np.linspace(25,end_P_1_3K,100)

A_3K_1,B_3K_1,C_3K_1 = curve_fit(expdec,P_3000_1,D_3000_1,p0=[0.01,0.6,30])[0]

plt.plot(P_3000_1,D_3000_1,"r-o")
plt.plot(P_fic1,expdec(P_fic1,A_3K_1,B_3K_1,C_3K_1))
plt.plot(P_3000,D_3000,"g-o")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()

end_P_2_3K=49
P_fic2_3K=np.linspace(end_P_2_3K,end_P_2,100)


P_3000_2=P_3000[P_3000>end_P_2_3K]
D_3000_2=D_3000[P_3000>end_P_2_3K]


A_3K_2,B_3K_2,C_3K_2 = curve_fit(expdec,P_3000_2,D_3000_2,p0=[0.05,0.78,50])[0]

plt.plot(P_3000_1,D_3000_1,"ro")
plt.plot(P_3000_2,D_3000_2,"bo")
plt.plot(P_fic1_3K,expdec(P_fic1,A_3K_1,B_3K_1,C_3K_1),"r-")
plt.plot(P_fic2_3K,expdec(P_fic2,A_3K_2,B_3K_2,C_3K_2),"b-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()


#=========================================================================

P_2000_1=P_2000[P_2000<end_P_1]
D_2000_1=D_2000[P_2000<end_P_1]


A_2K_1,B_2K_1,C_2K_1 = curve_fit(expdec,P_2000_1,D_2000_1,p0=[0.01,0.3,30])[0]

plt.plot(P_2000_1,D_2000_1,"ro")
plt.plot(P_fic1,expdec(P_fic1,A_2K_1,B_2K_1,C_2K_1),"r-")
plt.plot(P_2000,D_2000,"g-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()

P_2000_2=P_2000[P_2000>end_P_1]
D_2000_2=D_2000[P_2000>end_P_1]


A_2K_2,B_2K_2,C_2K_2 = curve_fit(expdec,P_2000_2,D_2000_2,p0=[0.01,0.6,50])[0]

plt.plot(P_2000_1,D_2000_1,"ro")
plt.plot(P_2000_2,D_2000_2,"bo")
plt.plot(P_fic1,expdec(P_fic1,A_2K_1,B_2K_1,C_2K_1),"r-")
plt.plot(P_fic2,expdec(P_fic2,A_2K_2,B_2K_2,C_2K_2),"b-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()


plt.plot(P_fic1,expdec(P_fic1,A_3K_1,B_3K_1,C_3K_1),"r-")
plt.plot(P_fic2,expdec(P_fic2,A_3K_2,B_3K_2,C_3K_2),"r-")
plt.plot(P_fic1,expdec(P_fic1,A_2K_1,B_2K_1,C_2K_1),"b-")
plt.plot(P_fic2,expdec(P_fic2,A_2K_2,B_2K_2,C_2K_2),"b-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()

#===========================================================================

P_2500_1=P_2500[P_2500<end_P_1]
D_2500_1=D_2500[P_2500<end_P_1]


A_25K_1,B_25K_1,C_25K_1 = curve_fit(expdec,P_2500_1,D_2500_1,p0=[0.01,0.3,30])[0]

plt.plot(P_2500_1,D_2500_1,"ro")
plt.plot(P_fic1,expdec(P_fic1,A_25K_1,B_25K_1,C_25K_1),"r-")
plt.plot(P_2500,D_2500,"g-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()

P_2500_2=P_2500[P_2500>end_P_1]
D_2500_2=D_2500[P_2500>end_P_1]


A_25K_2,B_25K_2,C_25K_2 = curve_fit(expdec,P_2500_2,D_2500_2,p0=[0.01,0.2,40])[0]

plt.plot(P_2500_1,D_2500_1,"ro")
plt.plot(P_2500_2,D_2500_2,"bo")
plt.plot(P_fic1,expdec(P_fic1,A_25K_1,B_25K_1,C_25K_1),"r-")
plt.plot(P_fic2,expdec(P_fic2,A_25K_2,B_25K_2,C_25K_2),"b-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.show()


plt.plot(P_3000,D_3000*5e-5,"ro")
plt.plot(P_fic1,expdec(P_fic1,A_3K_1,B_3K_1,C_3K_1)*5e-5,"r-")
plt.plot(P_fic2,expdec(P_fic2,A_3K_2,B_3K_2,C_3K_2)*5e-5,"r-")
plt.plot(P_2500,D_2500*5e-5,"go")
plt.plot(P_fic1,expdec(P_fic1,A_25K_1,B_25K_1,C_25K_1)*5e-5,"g-")
plt.plot(P_fic2,expdec(P_fic2,A_25K_2,B_25K_2,C_25K_2)*5e-5,"g-")
plt.plot(P_2000,D_2000*5e-5,"bo")
plt.plot(P_fic1,expdec(P_fic1,A_2K_1,B_2K_1,C_2K_1)*5e-5,"b-")
plt.plot(P_fic2,expdec(P_fic2,A_2K_2,B_2K_2,C_2K_2)*5e-5,"b-")
plt.xlabel("P (GPa)")
plt.ylabel("D (cm^2/s)")
plt.legend(["3000K","","","2500K","","","2000K","",""])
plt.show()


N=200
P_range=np.linspace(25,70,N)
file = open(str(folder+"Diffusion.dat"),"w") 
for i in range(N):
    if P_range[i] < end_P_1:
        file.write("2000 "+str(P_range[i])+" "+str(expdec(P_range[i],A_2K_1,B_2K_1,C_2K_1))+"\n")
    else:
        file.write("2000 "+str(P_range[i])+" "+str(expdec(P_range[i],A_2K_2,B_2K_2,C_2K_2))+"\n")
file.write("\n")
for i in range(N):
    if P_range[i] < end_P_1:
        file.write("2500 "+str(P_range[i])+" "+str(expdec(P_range[i],A_25K_1,B_25K_1,C_25K_1))+"\n")
    else:
        file.write("2500 "+str(P_range[i])+" "+str(expdec(P_range[i],A_25K_2,B_25K_2,C_25K_2))+"\n")
file.write("\n")     
for i in range(N):
    if P_range[i] < end_P_1:
        file.write("3000 "+str(P_range[i])+" "+str(expdec(P_range[i],A_3K_1,B_3K_1,C_3K_1))+"\n")
    else:
        file.write("3000 "+str(P_range[i])+" "+str(expdec(P_range[i],A_3K_2,B_3K_2,C_3K_2))+"\n")
file.write("\n")
file.close()


