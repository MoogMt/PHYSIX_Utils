#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 04:29:13 2020

@author: mathieumoog
"""


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
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump

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
            dist = getDistanceOrtho( positions, atom, atom2, cell_lengths )
            matrix[atom,atom2] = dist
            matrix[atom2,atom] = dist
    return matrix
def computeTransitionMatrix( states, nb_states, tau, step_max ):
    nb_step = len(states)
    matrix = np.zeros((nb_states,nb_states))
    for step in range( nb_step-step_max ):
        matrix[ states[step], states[step+tau] ] += 1
    return matrix
def computeChapmanKolmogorov( matrix, nb_states ):
    matrix_ck = np.zeros((nb_states,nb_states),dtype=float)
    for state_i in range( nb_states ):
        for state_j in range( nb_states ):
            for state_k in range( nb_states) :
                matrix_ck[ state_i, state_j ] +=  matrix[state_i,state_k]*matrix[state_k,state_j] 
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
max_neigh=5
nb_step=len(traj[:,0,0])

distances = np.zeros( (nb_step,nbC,max_neigh), dtype=float )
for step in range(nb_step):
    matrix = computeDistanceMatrix( traj[step,:,:], cell_lengths )
    distances[step,0:nbC,0:max_neigh] = np.sort(matrix,axis=1)[0:nbC,1:max_neigh+1]
  
    
nbins=25
r_width=0.5
plt.figure()      
plt.hist( np.reshape(distances[:,:,0],(20000*32)), bins=nbins, rwidth=r_width )
plt.hist( np.reshape(distances[:,:,1],(20000*32)), bins=nbins, rwidth=r_width )
plt.hist( np.reshape(distances[:,:,2],(20000*32)), bins=nbins, rwidth=r_width )
plt.hist( np.reshape(distances[:,:,3],(20000*32)), bins=nbins, rwidth=r_width )  
plt.show()    

nb_box=5
min_distance = np.zeros(max_neigh, dtype=float)
max_distance = np.zeros(max_neigh, dtype=float)
delta_box=np.zeros( max_neigh, dtype=float )
for neigh in range(max_neigh):
    min_distance[neigh] = np.min( distances[:,:,neigh] )
    max_distance[neigh] = np.max( distances[:,:,neigh] )
    delta_box[neigh] = (max_distance[neigh]-min_distance[neigh])/nb_box
    
    
states = np.zeros( (nb_step,nbC,max_neigh), dtype=int )
for neigh in range(max_neigh):
    states[:,:,neigh] = (distances[:,:,neigh]-min_distance[neigh])/delta_box[neigh]


msm = MarkovStateModel( lag_time=1, n_timescales=6)
msm.fit(states[:,0,0])
msm.timescales_


dt=5*0.001
frac = 0.5
max_step=int(nb_step*frac)
nb_tau_min=int(500)
nb_tau_max=int(2*nb_tau_min)

target_neigh=3

# Computing Transition Matrix for a given tau
matrix_tot=np.zeros((nb_box,nb_box,nb_tau_max), dtype=float )
matrix_tot_ck=np.zeros((nb_box,nb_box,nb_tau_min), dtype=float )
for tau in range(nb_tau_max):
    matrix = np.zeros((nb_box,nb_box),dtype=float)
    for carbon in range(nbC):
        matrix += computeTransitionMatrix( states[:,carbon,target_neigh-1], nb_box, tau, max_step )
    for state in range(nb_box):
        matrix[state,:] /= sum( matrix[state,:] )
    matrix_tot[:,:,tau] = matrix[:,:]
    if tau < nb_tau_min:
        matrix_tot_ck[:,:,tau] = computeChapmanKolmogorov( matrix_tot[:,:,tau], nb_box)  

plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("P_ij, P_ij^CK")
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[0,0,:], "k-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[0,0,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[1,1,:], "r-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[1,1,:], "r--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[2,2,:], "b-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[2,2,:], "b--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[3,3,:], "g-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[3,3,:], "g--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[4,4,:], "m-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[4,4,:], "m--" )
plt.show()    

target_state=9
plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("P_ij, P_ij^CK")
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[target_state,target_state,:], "k-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[target_state,target_state,:], "k--" )
plt.show()  

