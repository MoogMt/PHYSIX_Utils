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
max_neigh=5
nb_step=len(traj[:,0,0])
cut_off = 1.75
min_stat=1000

# Build States
coordC = np.zeros( (nb_step,nbC), dtype=int )
coordO = np.zeros( (nb_step,nbO), dtype=int )
for step in range(nb_step):
    matrix = computeContactMatrix( traj[step,:,:], cell_lengths, cut_off)
    for carbon in range(0,nbC):
        coordC[ step, carbon ] = int( sum(matrix[carbon,:]) )
    for oxygen in range(nbC,nb_atoms):
        coordO[ step, oxygen-nbC ] = int( sum(matrix[oxygen,:]) )
        
c_min = coordC.min()
o_min = coordO.min()
# Adapting the labels to make sure they are in the 0-nb_states range
coordC -= c_min
coordO -= c_min

msm = MarkovStateModel( lag_time=1, n_timescales=6)
msm.fit( coordC[:,0] )
msm.timescales_

# Computing nb of states (max)
nb_states_C = coordC.max()+1
nb_states_O = coordO.max()+1

# Computing Equilibrium States Probabilities
coordC_hist = np.zeros( nb_states_C )
ones_ = np.ones((nb_step,nbC), dtype=int )
for i in range( nb_states_C ):
    coordC_hist[i] = sum( ones_[ coordC == i ] )
# Clean marginal states
# for state in range( nb_states_C ):
#     if coordC_hist[state] < min_stat:
#         mask_to_clean = coordC[ :, : ]
coordC_hist /= sum(coordC_hist[:])

# Computing Equilibrium States Probabilities, cleaning marginals
ones_ = np.ones((nb_step,nbO), dtype=int )
coordO_hist = np.zeros( nb_states_O )
for i in range( nb_states_O ):
    coordO_hist[i] = sum( ones_[ coordO == i ] )
coordO_hist /= sum(coordO_hist[:])


# Plotting Oxygens
plt.figure()
plt.plot(coordC_hist,"b.-")
plt.plot(coordO_hist,"r.-")
plt.legend(["C states","O states"])
plt.show() 

dt=5*0.001
frac = 0.75
max_step=int(nb_step*frac)
nb_tau_min=int(250)
nb_tau_max=int(2*nb_tau_min)

# Computing Transition Matrix for a given tau 
matrix_tot=np.zeros((nb_states_C,nb_states_C,nb_tau_max), dtype=float )
matrix_tot_ck=np.zeros((nb_states_C,nb_states_C,nb_tau_min), dtype=float )
for tau in range(nb_tau_max):
    matrix = np.zeros((nb_states_C,nb_states_C),dtype=float)
    for carbon in range(nbC):
        matrix += computeTransitionMatrix( coordC[:,carbon], nb_states_C, tau+1, max_step )
    for state in range(nb_states_C):
        matrix[state,:] /= sum( matrix[state,:] )
    matrix_tot[:,:,tau] = matrix[:,:]
    if tau < nb_tau_min:
        matrix_tot_ck[:,:,tau] = computeChapmanKolmogorov( matrix_tot[:,:,tau], nb_states_C )  


carbon_target=3
matrix_markov = np.zeros( (4,4,nb_tau_min), dtype=float )
matrix_markov_ck = np.zeros( (4,4,nb_tau_min), dtype=float )
for tau in range(1,nb_tau_min+1):
    msm_matrix = MarkovStateModel( lag_time=tau, reversible_type="mle" ,n_timescales=nb_states_C, ergodic_cutoff="on", sliding_window=True, verbose=True)
    msm_matrix.fit( coordC[:,carbon_target] )
    matrix_markov[:,:,tau-1] = msm_matrix.transmat_
    for state_i in range( len(matrix_markov) ):
        for state_j in range( len(matrix_markov) ):
            for i in range( len(matrix_markov) ):
                matrix_markov_ck[ state_i, state_j, tau-1 ] +=  matrix_markov[state_i,i,tau-1]*matrix_markov[i,state_j,tau-1]

# PyEMMA
lags = [1,5,10,15,20,50,100,200]
implied_timescales = pe.msm.its(dtrajs=coordC[:,carbon_target].tolist(),lags=lags)
pe.plots.plot_implied_timescales(implied_timescales,units='time-steps', ylog=False)

M = pe.msm.estimate_markov_model(dtrajs=coordC[:,carbon_target].tolist(), lag = 10 )
cktest = M.cktest(nsets=3)
cktplt = pe.plots.plot_cktest(cktest)


plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("P_ij, P_ij^CK")
# plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[0,0,:], "k-" )
# plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[0,0,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_min,dt*1), matrix_markov[0,0,:], "k-" )
plt.plot( np.arange(0,2*dt*nb_tau_min,dt*2), matrix_markov_ck[0,0,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[1,1,:], "r-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[1,1,:], "r--" )
plt.plot( np.arange(0,dt*nb_tau_min,dt*1), matrix_markov[0,1,:], "k-" )
plt.plot( np.arange(0,2*dt*nb_tau_min,dt*2), matrix_markov_ck[0,1,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[1,2,:], "b-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[1,2,:], "b--" )
plt.plot( np.arange(0,dt*nb_tau_min,dt*1), matrix_markov[0,2,:], "k-" )
plt.plot( np.arange(0,2*dt*nb_tau_min,dt*2), matrix_markov_ck[0,2,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[1,3,:], "g-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[1,3,:], "g--" )
plt.plot( np.arange(0,dt*nb_tau_min,dt*1), matrix_markov[0,3,:], "k-" )
plt.plot( np.arange(0,2*dt*nb_tau_min,dt*2), matrix_markov_ck[0,3,:], "k--" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*1), matrix_tot[1,4,:], "m-" )
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), matrix_tot_ck[1,4,:], "m--" )
plt.show()    
    
rmseC = np.zeros(nb_tau_min, dtype=float)
for tau in range(nb_tau_min):
    mat = matrix_tot[:,:,2*tau]-matrix_tot_ck[:,:,tau]
    rmseC[tau] = sum(sum( mat*mat ))/(nb_states_C*nb_states_C)
   
plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("RMSE C (%)")
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), rmseC*100 )
plt.show()

matrix_tot=np.zeros((nb_states_O,nb_states_O,nb_tau_max), dtype=float )
matrix_tot_ck=np.zeros((nb_states_O,nb_states_O,nb_tau_min), dtype=float )
for tau in range(nb_tau_max):
    matrix = np.zeros((nb_states_O,nb_states_O),dtype=float)
    for carbon in range(nbC):
        matrix += computeTransitionMatrix( coordO[:,carbon], nb_states_O, tau, max_step )
    for state in range(nb_states_O):
        matrix[state,:] /= sum( matrix[state,:] )
    matrix_tot[:,:,tau] = matrix[:,:]
    if tau < nb_tau_min:
        matrix_tot_ck[:,:,tau] = computeChapmanKolmogorov( matrix_tot[:,:,tau], nb_states_O )  


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
plt.show()    

rmseO = np.zeros(nb_tau_min, dtype=float)
for tau in range(nb_tau_min):
    mat = matrix_tot[:,:,2*tau]-matrix_tot_ck[:,:,tau]
    rmseO[tau] = sum(sum( mat*mat ))/(nb_states_O*nb_states_O)
   
plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("RMSE O (%)")
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), rmseO*100 )
plt.show()

plt.figure()
plt.xlabel("Time lag (ps)")
plt.ylabel("RMSE all (%)")
plt.plot( np.arange(0,dt*nb_tau_max,dt*2), (rmseO+rmseC)*100*0.5 )
plt.show()


