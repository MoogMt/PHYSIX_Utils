#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Fonctions Radiales

#------------------------------------------------

#------------------------------------------------
"""

#-------------------
# Bibliotheque plot
#---------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
#---------------------------------------------

# Bibliotheques de calcul
#---------------------------------------------
import numpy as np
# Factorielle
from scipy.special import factorial
#---------------------------------------------

#----------------------
# Polynome de Laguerre
#---------------------------------------------
def laguerre(n,p,x) :
    if n == 0 :
        return 1.
    elif n==1 :
        return 1. + float(p) - float(x)
    else :
        lm = 1.;
        l = 1. + float(p) - float(x)
        for k in range(2,n+1):
            lp = ( (2.*k + 1. + p - x )*l - (k+p)*lm )/( k + 1.) ;
            lm = float( l ) ;
            l = float( lp ) ;
        return float(l)
#---------------------------------------------

#--------------
# Start Message
#---------------------------------------
print(" ")
print("---------------------------------");
print("| Calcul des fonctions radiales |")
print("---------------------------------");
print(" ")
#---------------------------------------

#--------------------------------------
# Getting the principal quantum number
#--------------------------------------
n = input("Enter n: ");
#--------------------------------------

#-------------------
# Physical variables
#-------------------------------
a0 = 1.00# Bohr Radius
#-------------------------------

#------
# Step
#-------------------------------
r_step = 300
r_max = 10.0
dr = r_max/float(r_step)
#-------------------------------

#-------------------------------
print(" ")
print("Calcul des fonctions radiales...")
print(" ")
#-------------------------------

#-------------
# Compute c(l)
#----------------------------------------------------------------------------------------------
c = np.zeros(( n ))
for l in range( n ):
    c[l] = float(factorial( n-l-1 ))/( 2*n* float(factorial(n+l))**3 )*( 2/(n*a0) )**(l+1.5)
#----------------------------------------------------------------------------------------------

#------------------------------------
# Boucle sur toutes les valeurs de r
#------------------------------------------------------------------------------------------
R = np.zeros(( n , r_step+1 ))
for i in range( r_step+1) :
    # Compute r
    r = i*dr
    # Boucle qui couvre toutes les valeurs de l (0 < l < n)
    for l in range( n ):
        # Calcul de la fonction radial pour r
        R[l,i] = c[l] * r**l * np.exp( -r/(n*a0) ) * laguerre( n-l-1 , 2*l+1 , 2*r/(n*a0) ) ;
#------------------------------------------------------------------------------------------

#----------------
# Normalisation
#------------------------------
norm = 0.0000000000
for l in range( n ):
    for i in range( r_step+1 ):
        norm += (i*dr*R[l,i])**2*dr
    R[l,:] /= np.sqrt(norm)
#------------------------------

#--------------------------------------------------
# Plot toutes les fonctions radiales pour n fixe
#--------------------------------------------------
r = np.linspace(0,r_max,r_step+1)
for l in range(n):
      plt.plot( r , R[l,:] , label="check")
      plt.xlabel('r (Bohr)')
      plt.ylabel('R_n(r)')
plt.show()
#---------------------------------------------------
      
#-------------------------------
print(" ")
print("Travail termine!")
print(" ")
#-------------------------------
