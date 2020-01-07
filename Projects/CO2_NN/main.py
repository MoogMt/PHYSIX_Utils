#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:50:29 2019

@author: moogmt
"""

import nnmetadata as mtd
import filexyz as xyz
import cpmd 
import descriptors as desc
import pandas as pd
import behler
import keras

data_base  = "/media/moogmt/Elements/CO2/"

verbose_check=True # Whether we want error messages 
debug = False # Debug mode: verbose descriptions are written in a debug file
# in order to check everything

volume=8.82
temperature=3000 
run_nb=1

folder_in = data_base + str(volume) + "/" + str(temperature) + "K/" + str(run_nb) + "-run/"
folder_out = data_base + str(volume) + "/" + str(temperature) + "K/Data/"

file_traj = folder_in + "TRAJEC.xyz"
file_energies = folder_in + "ENERGIES"


# EXTRACTING DATA
#=============================================================================#
metadata=mtd.buildMetaData(file_traj,file_energies,folder_out, temperature)
if not mtd.checkMetaDataIO(metadata,verbose_check):
    exit
nb_step=cpmd.getNbLineEnergies(file_energies)
# Reading trajectory
traj = xyz.readPbcCubic( file_traj, volume )
metadata['periodic'] = True
# Reading ENERGIES file
energies=cpmd.readPotEnergy( file_energies )
# Homogenizing with stride
stride_energies=5 # Improve by reading the stride in the input file and making adjustements.
energies=energies[0:len(energies):stride_energies]
# add a check to verify congruence of sizes...
# Getting species present in the simulation
metadata['n_atoms'] = len(traj[0])
metadata['species'] = mtd.getSpecies(traj[0])
metadata['n_species'] = len(metadata['species'])
#traj=mtd.sortAtomsSpecie(traj) # Use if you need to sort the atoms per type, current implementation is *very* slow
metadata["species_sorted"]=True
metadata=mtd.getStartSpecies(traj,metadata)
metadata=mtd.getNbAtomsPerSpecies(traj,metadata)
#=============================================================================#

# CREATING DESCRIPTORS
#=============================================================================#
# Creating training set
metadata['n_jobs'] = 8 # Number of parallel cores to use (CPU)
metadata['train_set_size'] = 400
metadata['total_size_set'] = len(energies)
metadata,data_train = mtd.choseDataRandom(metadata,traj,energies)
# Creating testing set
metadata["test_set_size"] = 500
data_test = mtd.ChooseDataExclusion(metadata,traj,energies)
# Build descriptors from positions (train set only)
sigma_  = 0.9  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 3.5 # cut_off SOAP, 
nmax_   = 2 
lmax_   = 2
metadata, descriptors = desc.createDescriptorsSOAP(data_train,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
descriptors, scalers = mtd.scaleData(descriptors,metadata)
data_train=data_train.join(pd.DataFrame({'descriptor':list(descriptors)}))
metadata, descriptors = desc.createDescriptorsSOAP(data_test,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
data_test=data_test.join(pd.DataFrame({'descriptor':list(descriptors)}))
#=============================================================================#

# BUILDING NETWORK
#=============================================================================#
model=behler.buildNetwork(metadata)
model.compile(loss=metadata["loss_func"], optimizer=metadata["optimizer"], metrics=['accuracy'])
if metadata["plot_network"]:
    keras.utils.plot_model(model,to_file=metadata["path_plot_network"])
#=============================================================================#

# TRAINING NETWORK
#=============================================================================#
data_train,data_test,mean_error,metadata = behler.train(data_train,data_test,model,metadata)
#=============================================================================#

