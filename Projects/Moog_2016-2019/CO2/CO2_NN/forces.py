
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

path_sim = str( "/Users/mathieumoog/Documents/CO2/" + str(volume) + "/" + str(temperature) + "K/" )

cell_lengths = np.ones(3)*volume

traj_path = str( path_sim + "TRAJEC_fdb_wrapped.xyz" )
traj = filexyz.readAsArray( traj_path )

nbC=32
nbO=64
nb_atoms=nbC+nbO
max_neigh=4
nb_step=len(traj[:,0,0])
cut_off = 1.75