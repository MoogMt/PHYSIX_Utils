#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 21:31:51 2020

@author: moogmt

The goal here is to provide an in-depth analysis of CO2 under extreme conditions 
using a combination of:
    - Very potent local descriptor with the Smooth Overlap of Atomic Positions (SOAP)
    - A Feed-Forward Neural Network 
    - Markov State Models
    
Here are dealing with a highly reactive polymeric fluid where bonds are not easily
described. 

The code functions as follows:
    - We first label all carbons depending on their bonding coordination:
    carbon that have neighbors in the contentious region (1.6 < d < 2) are 
    excluded and the remaining carbons are labeled depending on the number of
    neighbors that they have (C/O separately)
    - We then construct simple first and second layer atomic states
    - The first two steps are repeated with oxygen
    - SOAP descriptors are then computed for all atoms, and the resulting
    vectors are separated into a training and test set and fed to a NN that
    aims at giving them the proper label. Incase some labels are too weak, they 
    are bootstrapped so that there is no distribution issue.
    - The same procedure is followed with a SVM machine for classification and
    results are compared
    - The best procedure is then used to label the C/O that could not be labelled
    earlier. A handshake is used to make sure that the labellign give consistent 
    results. If not a more clever method should be used.
    - Second layers states are then computed and and used as reference for a 
    Markov State Model.
    - If the results are conclusive, the lifetime of states are computed
    for each states and if possible the mean first passage time from each structure
    to each other as well.
    - Each state is characterized as a chemical active site and the lengths and angles
    dsitrbutions are computed.
"""

import nnmetadata as mtd
import filexyz as xyz
import cpmd 
import descriptors as desc
import numpy as np

data_base  = "/media/moogmt/Elements/CO2/"

verbose_check=True # Whether we want error messages 
debug = False # Debug mode: verbose descriptions are written in a debug file
# in order to check everything


# EXTRACTING DATA
#=============================================================================#
volume=8.82
temperature=3000 
run_nb=1
folder_in = data_base + str(volume) + "/" + str(temperature) + "K/" + str(run_nb) + "-run/"
folder_out = data_base + str(volume) + "/" + str(temperature) + "K/Data/"
file_traj = folder_in + "TRAJEC.xyz"
file_energies = folder_in + "ENERGIES"
#------------------------------------------------------------------------------
metadata=mtd.buildMetaData(file_traj,file_energies,folder_out, temperature)
if not mtd.checkMetaDataIO(metadata,verbose_check):
    exit
#------------------------------------------------------------------------------
nb_step=cpmd.getNbLineEnergies(file_energies)
# Reading trajectory
traj = xyz.readPbcCubic( file_traj, volume )
periodic = True
# Reading ENERGIES file
data_out     = cpmd.readEnergiesFile( file_energies )
energies     = cpmd.extractPotentialEnergy(     data_out )
comput_time  = cpmd.extractSCFcomputationTime(  data_out )
# Homogenizing with stride
stride_energies=5 # Improve by reading the stride in the input file and making adjustements.
energies    = energies    [ 0: len(energies):stride_energies ]
comput_time = comput_time [ 0: len(energies):stride_energies ]
# add a check to verify congruence of sizes...
# Getting species present in the simulation
n_atoms = len(traj[0])
species = mtd.getSpecies(traj[0])
n_species = len(species)
#traj=mtd.sortAtomsSpecie(traj) # Use if you need to sort the atoms per type, current implementation is *very* slow
species_sorted=True
total_size_set = len( energies )
start_species=mtd.getStartSpecies( traj, species )
nb_element_species=mtd.getNbAtomsPerSpecies( traj, species )
#=============================================================================#

# Keeping only easily identifiable carbons
#=============================================================================#
import ase.geometry as asegeom
import matplotlib.pyplot as plt

# Build descriptors from positions (train set only)
sigma_  = 0.9  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 4.0 # cut_off SOAP, 
nmax_   = 3
lmax_   = 2

specie_labeled=np.zeros((0,2))
specie=0
step=0
distances=asegeom.get_distances(traj[step].positions, pbc=traj[step].pbc,cell=traj[step].cell )[1]
label=[ sum() for atom in range(len(traj)) ]
index_=[ sum( dist > 1.6 and dist < 1.8 for dist in distances[atom,:] ) < 1 for atom in range(len(traj[0])) ]
descriptorC=desc.createDescriptorsSOAP( traj[step], species, sigma_, cutoff_, nmax_, lmax_, periodic )[0:32][index_[0:32],:]
descriptorO=desc.createDescriptorsSOAP( traj[step], species, sigma_, cutoff_, nmax_, lmax_, periodic )[0:32][index_[0:32],:]


#=============================================================================#

