# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#===================================
import numpy as np
import scipy as scp
from scipy.misc import factorial
import matplotlib.pyplot as plt
#===================================

#==================
# Fonctions utiles
#======================================================================
# Calcul la fonction radial
# Lien: 
def laguerre(n,p,x):
    if n == 0:
        return 1.0
    elif n == 1:
        return 1+p-x
    else:
        lm = 1
        l =1+p-x
        for k in range(2,n+1):
            lp = ((2*k+1+p-x)*l-(k+p)*lm)/(k+1)
            lm = l
            l=lp
        return l
    return laguerre
#---------------------------------------------------------------------
# Calcule la fonction radiale
# Lien:
def radial(r,n,l):
    # Coefficient c(l)
    c_l= factorial(n-l-1)/(2*n*factorial(n+1)**3)*(2/(n*a0))**(l+3/2.)
    # Retourne l'expression complète
    return c_l*r**l*np.exp(-r/(n*a0))*laguerre(n-l-1,2*l+1,2*r/(n*a0))
def radial_trace(r,n,l,dr):
    # On initialise le vecteur pour la fonction radiale
    radial_=np.ones(len(r))
    # Facteur de normalisation
    norm=0
    # Boucle de r=0 à r=rmax
    for i in range(len(r)):
        # Calcul de la fonction radiale
        radial_[i]=radial(i*dr,n,l);
        # On rajoute la valeur au facteur de normalisation
        norm+=abs(radial_[i]);
    # On renvoit le résultat
    return radial_/norm
#======================================================================

#=============
# Constantes
#============================
# Rayon de Bohr, en Bohr
a0=1.000000
#=============================

# Valeurs pour plot
#===================================================================
rmax=(float)(input("Donner la valeur max de r à plotter (Bohr): "))
nb_r=(int)(input("Donner le nombre de points: "))
dr=rmax/float(nb_r)
#===================================================================

# Calcul de la fonction radiale
#=====================================
r_list=np.linspace(0,rmax,nb_r);
#=====================================

# Plots pour plusieurs n
#==============================================
plt.figure(1)
plt.plot(r_list,radial_trace(r_list,1,0,dr))
plt.plot(r_list,radial_trace(r_list,2,0,dr))
plt.plot(r_list,radial_trace(r_list,3,0,dr))
plt.plot(r_list,radial_trace(r_list,4,0,dr))
plt.xlabel("r (Bohr)")
plt.ylabel("R_n,l(r)")
plt.legend(["n=1,l=0","n=2, l=0","n=3, l=0","n=4,l=0"])
plt.show()
#=============================================

# Plots pour n fix et tout les l
#=============================================
plt.figure(2)
plt.plot(r_list,radial_trace(r_list,3,0,dr))
plt.plot(r_list,radial_trace(r_list,3,1,dr))
plt.plot(r_list,radial_trace(r_list,3,2,dr))
plt.xlabel("r (Bohr)")
plt.ylabel("R_n,l(r)")
plt.legend(["n=3,l=0","n=3, l=1","n=3, l=2"])
plt.show()
#=============================================
