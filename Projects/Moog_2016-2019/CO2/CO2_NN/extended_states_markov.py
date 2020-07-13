#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 05:54:11 2020

@author: mathieumoog
"""

import cpmd
import filexyz
import numpy as np
import matplotlib.pyplot as plt
# MSMbuilder ( lacks CK validation )
from msmbuilder.msm import MarkovStateModel
from msmbuilder.msm import BayesianMarkovStateModel
from msmbuilder.utils import dump
# PyEMMMA ( has CK validation )
import pyemma as pe
from pyemma.datasets import double_well_discrete

def getDistance1Dsq( position1, position2, length):
    dist = position1-position2
    half_length = length*0.5
    if dist > half_length :
        dist -= length
    elif dist < -half_length:
        dist += length
    return dist*dist
def getDistanceOrtho( positions, index1, index2, cell_lengths ):
    dist=0
    for i in range(3):
        dist += getDistance1Dsq( positions[index1,i], positions[index2,i], cell_lengths[i] )
    return np.sqrt(dist)
def computeContactMatrix( positions, cell_lengths, cut_off ):
    nb_atoms = len(positions[:,0])
    matrix = np.zeros(( nb_atoms, nb_atoms ))
    for atom in range(nb_atoms):
        for atom2 in range(atom+1,nb_atoms):
            if getDistanceOrtho( positions, atom, atom2, cell_lengths ) < cut_off :
                matrix[atom,atom2] = 1
                matrix[atom2,atom] = 1 
    return matrix
def computeDistanceMatrix( positions, cell_lengths ):
    nb_atoms = len(positions[:,0])
    matrix = np.zeros(( nb_atoms, nb_atoms ))
    for atom in range(nb_atoms):
        for atom2 in range(atom+1,nb_atoms):
            matrix[atom,atom2] = getDistanceOrtho( positions, atom, atom2, cell_lengths )
            matrix[atom2,atom] = matrix[atom,atom2]
    return matrix
def computeTransitionMatrix( states, nb_states, tau, step_max ):
    nb_step = states.shape[0]
    matrix = np.zeros((nb_states,nb_states))
    for step in range( nb_step-step_max ):
        matrix[ int(states[step]), int(states[step+tau]) ] += 1
    return matrix
def computeTransitionMatrixClean( states, nb_states, tau, step_max ):
    nb_step = states.shape[0]
    matrix = np.zeros((nb_states,nb_states))
    for step in range( nb_step-step_max ):
        start = int( states[step] )
        end   = int( states[step+tau] )
        if start >= 0 and end >= 0:
            matrix[ start, end ] += 1
    return matrix
def computeChapmanKolmogorov( matrix, nb_states ):
    matrix_ck = np.zeros((nb_states,nb_states),dtype=float)
    for state_i in range( nb_states ):
        for state_j in range( nb_states ):
            for i in range(nb_states):
                matrix_ck[ state_i, state_j ] += matrix[state_i,i]*matrix[i,state_j] 
    return matrix_ck



volume=8.82
temperature=3000
# run_nb=1

path_sim = str( "/Users/mathieumoog/Documents/CO2/" + 
               str(volume) + "/" + 
               str(temperature) + "K/" 
               # + str(run_nb) + "-run/" 
               )


cell_lengths = np.ones(3)*volume

traj_path = str( path_sim + "TRAJEC_fdb_wrapped.xyz" )
traj = filexyz.readAsArray( traj_path )

nbC=32
nbO=64
nb_atoms=nbC+nbO
max_neigh=4
nb_step=len(traj[:,0,0])
cut_off = 1.75


states=[]
states_time = np.ones( ( nb_step, nbC ), dtype=float )*(-1)
# Loop over states 
for step in range(nb_step):
    # Compute all distances between all atoms
    matrix = computeDistanceMatrix( traj[step,:,:], cell_lengths )
    # Loop over carbon atoms
    for carbon in range(nbC):
        local_state = np.ones( max_neigh )*(-1)
        # Get the index of the nearest oxygen neighbors of carbon target
        neigh_index = np.argsort( matrix[ carbon, nbC:nb_atoms ] )[0:max_neigh]
        # Loop over each of the neighbor of the target carbon
        for neigh in range(max_neigh):
            # If the neighbor is linked with the target carbon 
            if matrix[ carbon, neigh_index[neigh]+nbC ] < cut_off :
                # Counting the number of neighbors that each neighbor has (remove 1 for self)
                local_state[neigh] = sum( np.ones(nb_atoms)[ matrix[ neigh_index[neigh]+nbC, : ] < cut_off]  ) - 1
        # Allocate the state number
        for state_index, state in enumerate(states):
            if sum( state -  local_state ) == 0 :
                states_time[ step, carbon ] = state_index
        # If the state is not found we add the state to the database
        if states_time[ step, carbon ] == -1 : 
            states.append( local_state )
            states_time[ step, carbon ] = len( states ) - 1
states = np.array( states )    
         
nb_states=len(states)
   
ones_ = np.ones((nb_step,nbC), dtype=int )
states_hist = np.zeros( nb_states )
for state in range( nb_states ):
    states_hist[state] = sum( ones_[ states_time == state ] )
states_hist /= sum( states_hist[:] )


# Plotting Oxygens
plt.figure()
plt.plot(states_hist,"b.-")
plt.legend(["C states"])
plt.show()      
          
dt=5*0.001
frac = 0.75
max_step=int(nb_step*frac)
nb_tau_min=int(250)
nb_tau_max=int(2*nb_tau_min)

# Computing Transition Matrix for a given tau 
matrix_tot=np.zeros((nb_states,nb_states,nb_tau_max), dtype=float )
matrix_tot_ck=np.zeros((nb_states,nb_states,nb_tau_min), dtype=float )
for tau in range(nb_tau_max):
    matrix = np.zeros((nb_states,nb_states),dtype=float)
    for carbon in range(nbC):
        matrix += computeTransitionMatrix( states_time[:,carbon], nb_states, tau+1, max_step )
    for state in range(nb_states):
        matrix[state,:] /= sum( matrix[state,:] )
    matrix_tot[:,:,tau] = matrix[:,:]
    if tau < nb_tau_min:
        matrix_tot_ck[:,:,tau] = computeChapmanKolmogorov( matrix_tot[:,:,tau], nb_states )    
            
plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("P_ij, P_ij^CK")
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[0,0,:],    "r-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[0,0,:], "r--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[1,1,:],    "b-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[1,1,:], "b--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[2,2,:],    "g-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[2,2,:], "g--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[3,3,:],    "m-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[3,3,:], "m--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[4,4,:],    "k-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[4,4,:], "k--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[5,5,:],    "c-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[5,5,:], "c--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[6,6,:],    "m-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[6,6,:], "m--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot[7,7,:],    "g-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_ck[7,7,:], "g--" )
plt.show()    
            

cut_off_percent = 0.025
states_time_cleaned = states_time
for step in range( nb_step ):
    for carbon in range(nbC):
        if states_hist[ int(states_time[step,carbon]) ] < cut_off_percent :
            states_time_cleaned[step,carbon] = - 1 

nb_states_cleaned=0
for state in range(nb_states):
    if states_hist[state] > cut_off_percent :
        nb_states_cleaned += 1

matrix_tot_clean    = np.zeros( ( nb_states_cleaned, nb_states_cleaned, nb_tau_max ), dtype=float )
matrix_tot_clean_ck = np.zeros( ( nb_states_cleaned, nb_states_cleaned, nb_tau_min ), dtype=float )
for tau in range(nb_tau_max):
    matrix = np.zeros( ( nb_states_cleaned, nb_states_cleaned ), dtype=float )
    for carbon in range(nbC):
        matrix += computeTransitionMatrixClean( states_time_cleaned[:,carbon], nb_states_cleaned, tau+1, max_step )
    for state in range(nb_states_cleaned):
        matrix[state,:] /= sum( matrix[state,:] )
    matrix_tot_clean[:,:,tau] = matrix[:,:]
    if tau < nb_tau_min:
        matrix_tot_clean_ck[:,:,tau] = computeChapmanKolmogorov( matrix_tot_clean[:,:,tau], nb_states_cleaned )    

plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("P_ij, P_ij^CK")
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[0,0,:],    "r-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[0,0,:], "r--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[1,1,:],    "b-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[1,1,:], "b--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[2,2,:],    "g-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[2,2,:], "g--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[3,3,:],    "m-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[3,3,:], "m--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[4,4,:],    "k-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[4,4,:], "k--" )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*1 ), matrix_tot_clean[5,5,:],    "c-"  )
plt.plot( np.arange( 0, dt*nb_tau_max, dt*2 ), matrix_tot_clean_ck[5,5,:], "c--" )
plt.show()   

#---------------------------------------------------------------------------------------------------

# Cleaning memory
matrix_tot=[]
matrix_tot_tot=[]
matrix_tot_ck=[]
matrix_tot_clean=[]
matrix_tot_clean_ck=[]
states_time=[]


