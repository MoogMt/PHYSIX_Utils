#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
# Harmoniques spheriques

# Ce code contient des informations permettant de tracer 
# des harmoniques sphériques.
"""

#-------------------
# Bibliotheque plot
#---------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
# Plot 3d
from mpl_toolkits.mplot3d import Axes3D
#---------------------------------------------

# Bibliotheques de calcul
#---------------------------------------------
import numpy as np
# Polynome de Legendre
#from scipy.special import legendre
# Polynome associe de Legendre
from scipy.special import lpmn
# Factorielle
from scipy.special import factorial
#---------------------------------------------

#-----------------------------------
# Calcul de l'harmonique spherique
#------------------------------------------------------------------------------
def calculHarmoniqueSpherique( l , m , theta, phi):
    return np.sqrt(2*factorial(l-m)/factorial(l+m))*lpmn(m,l,np.cos(theta))[0][m][l]*np.cos(m*phi)
#------------------------------------------------------------------------------

#------------------------------
# Ecriture du message de debut
#---------------------------------------------------
print(" ")
print("=====================================")
print("| Calcul des harmoniques spheriques |")
print("=====================================")
print(" ")
#---------------------------------------------------

#------------------------
# Recuperation de l et m
#------------------------------------
l = input("Entrez une valeur de l: ")
m = input("Entrez une valeur de m: ")
if l < m :
    print("Merci d'entrer des valeurs de l,m telles que l > m !")
    exit;
graph_lab="Ylm(theta,phi) l={:g} m={:g} ".format(l,m)
#------------------------------------

# Affichage du message de travail
#------------------------------------------
print(" ")
print("Calcul des harmoniques en cours...")
print(" ")
7#------------------------------------------

#-----------------------------------
# Computes the spherical harmonics
#-----------------------------------------------------------------------------
dummy = 0;
step = 5;
# Tableaux pour angles
theta_tab = np.arange(0,180+1,step); phi_tab = np.arange(0,360+1,step)
# Tableaux pour les harmoniques spheriques
x = np.zeros(( len(theta_tab)*len(phi_tab) ))
y = np.zeros(( len(theta_tab)*len(phi_tab) ))
z = np.zeros(( len(theta_tab)*len(phi_tab) ))
# Calcul sur toutes les valeurs de theta et phi
for theta in theta_tab:
    for phi in phi_tab:
        # Calcul des Ylm(theta,phi)
        mod =  calculHarmoniqueSpherique( l , m , theta*np.pi/180., phi*np.pi/360.)
        x[dummy] = mod*np.cos(theta)*np.sin(phi)
        y[dummy] = mod*np.sin(theta)*np.sin(phi)
        z[dummy] = mod*np.cos(phi)
        dummy += 1
#----------------------------------------------------------------------------

#---------------------------------------------------------------------
print(" ")
print("Travail termine, fermer le graphe pour arreter le programme.")
print(" ")
#----------------------------------------------------------------------

#-----------------
# Plotting graph
#----------------------------------------------------
fig = plt.figure()
ax = fig.gca( projection="3d" )
ax.scatter( x, y, z , label=graph_lab )
ax.legend( graph_lab )
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
# Show
plt.show()
#----------------------------------------------------

