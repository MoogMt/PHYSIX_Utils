#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:50:29 2019

@author: moogmt inspired by the work of julienh
"""

import nnmetadata as mtd
import filexyz as xyz
import cpmd 
import descriptors as desc
import pandas as pd
import keras
import numpy as np

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
metadata['train_set_size'] = 2000
metadata['total_size_set'] = len(energies)
metadata,data_train = mtd.choseTrainDataRandom(metadata,traj,energies)
# Creating testing set
metadata["test_set_size"] = 1000
metadata, data_test = mtd.choseTestDataRandomExclusion(metadata,traj,energies)
# Build descriptors from positions (train set only)
sigma_  = 0.8  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 3.2 # cut_off SOAP, 
nmax_   = 2 
lmax_   = 1
# Train set
#-----------------------------------------------------------------------------
metadata, descriptors = desc.createDescriptorsSOAP(data_train,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
descriptors, scalers = mtd.scaleData(descriptors,metadata)
data_train=data_train.join(pd.DataFrame({'descriptor':list(descriptors)}))
# Scaling
# Test set
#------------------------------------------------------------------------------
metadata, descriptors = desc.createDescriptorsSOAP(data_test,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
data_test=data_test.join(pd.DataFrame({'descriptor':list(descriptors)}))
# Building Inputs
input_train=mtd.buildInput(data_train)
output_train=data_train.energy.values
input_test=mtd.buildInput(data_test)
output_test=data_test.energy.values
#=============================================================================#

# BUILDING NETWORK
#=============================================================================#
# Parameters of the Neural net
import behler
        
# Iteration parameters
metadata["loss_fct"] = 'mean_squared_error' # Loss function in the NN
metadata["optimizer"] = 'adam'                    # Choice of optimizers for training of the NN weights 
metadata["n_epochs"] = 500                  # Number of epoch for optimization?
metadata["patience"] = 100                  # Patience for convergence
metadata["restore_weights"] = True
    
# Subnetorks structure
metadata["activation_fct"] = 'tanh'  # Activation function in the dense hidden layers
metadata["n_nodes_per_layer"] = 30           # Number of nodes per hidden layer
metadata["n_hidden_layer"] = 5                # Number of hidden layers
metadata["n_nodes_structure"]=np.ones((metadata["n_species"],metadata["n_hidden_layer"]),dtype=int)*metadata["n_nodes_per_layer"] # Structure of the NNs (overrides the two precedent ones)
        
# Dropout coefficients
metadata["dropout_coef"]=np.zeros((metadata["n_species"],metadata["n_hidden_layer"]+1)) # Dropout for faster convergence (can be desactivated) 
metadata["dropout_coef"][0,:]=0.2
metadata["dropout_coef"][1:,:]=0.5
        
# Plot network
metadata["plot_network"]=False
metadata["path_plot_network"]=str(folder_out+"plot_network.png")
metadata["saved_model"] = False
metadata["path_folder_save"]=str(folder_out)
metadata["suffix_write"]=str("train-"+metadata["train_set_size"] + "_" +
           "test-"+metadata["test_set_size"]        + "_" +
           "layers-"+metadata["n_hidden_layer"]     + "_" +
           "n_nodes-"+metadata["n_nodes_per_layer"] + "_" + 
           "nmaxSOAP-" + metadata["nmax_SOAP"] + "_" + 
           "lmaxSOAP-" + metadata["lmax_SOAP"] + "_" +
           "sigmaSOAP-" + metadata["sigma_SOAP"] + "_"+
           "cutoffSOAP-" + metadata["cutoff_SOAP"] 
           )

metadata, metadata_stat, predictions_train, prediction_test = behler.buildTrainPredictWrite(metadata,input_train,input_test,output_train,output_test)

