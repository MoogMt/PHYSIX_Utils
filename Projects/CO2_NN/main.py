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
import numpy as np
import behler
import keras

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
metadata['train_set_size'] = 1000
metadata['total_size_set'] = len(energies)
metadata, input_train_raw, output_train = mtd.choseTrainDataRandom(metadata,traj,energies)
# Creating testing set
metadata["test_set_size"] = 200
metadata, input_test_raw, output_test = mtd.choseTestDataRandomExclusion(metadata,traj,energies)
# Build descriptors from positions (train set only)
sigma_  = 0.9  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 3.2 # cut_off SOAP, 
nmax_   = 3
lmax_   = 3
# Train set
#-----------------------------------------------------------------------------
metadata, input_train = desc.createDescriptorsSOAP(input_train_raw,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
# Test set
#------------------------------------------------------------------------------
metadata, input_test = desc.createDescriptorsSOAP(input_test_raw,metadata,sigma_SOAP=sigma_,cutoff_SOAP=cutoff_,nmax_SOAP=nmax_,lmax_SOAP=lmax_)
#------------------------------------------------------------------------------
# Scaling 
#------------------------------------------------------------------------------
# Scaling Energy
metadata, output_train_scale = mtd.scaleEnergy( output_train, metadata )
metadata, output_test_scale  = mtd.scaleEnergy( output_test,  metadata )
# Scaling Input
scalers = mtd.createScaler( input_train, metadata ) # Create scaler on training set
input_train_scale = mtd.applyScale( scalers, input_train, metadata )
input_test_scale  = mtd.applyScale( scalers, input_test,  metadata )
# PCA
#------------------------------------------------------------------------------
#metadata["path_pcavar"]=str(folder_out+"pca_var.dat")
#var_pca = mtd.pcaVariance(input_train,metadata["path_pcavar"])
#=============================================================================#

# BUILDING NETWORK
#=============================================================================#
# Parameters of the Neural net

# Iteration parameters
metadata["loss_fct"] = 'mean_squared_error' # Loss function in the NN
metadata["optimizer"] = 'Adam'                    # Choice of optimizers for training of the NN weights 
metadata["n_epochs"] = 1000                  # Number of epoch for optimization?
metadata["patience"] = 100                  # Patience for convergence
metadata["restore_weights"] = True
metadata["batch_size"] = 500
metadata["verbose_train"] = 1
metadata["early_stop_metric"]=['mse']

# Subnetorks structure
metadata["activation_fct"] = 'tanh'  # Activation function in the dense hidden layers
metadata["n_nodes_per_layer"] = 30           # Number of nodes per hidden layer
metadata["n_hidden_layer"] = 3               # Number of hidden layers
metadata["n_nodes_structure"]=np.ones((metadata["n_species"],metadata["n_hidden_layer"]),dtype=int)*metadata["n_nodes_per_layer"] # Structure of the NNs (overrides the two precedent ones)
metadata["kernel_constraint"] = keras.constraints.maxnorm(2)

# Dropout coefficients
metadata["dropout_coef"]=np.zeros((metadata["n_species"],metadata["n_hidden_layer"]+1)) # Dropout for faster convergence (can be desactivated) 
metadata["dropout_coef"][0,:]=0.2
metadata["dropout_coef"][1:,:]=0.5
    
# Plot network
metadata["plot_network"]=True
metadata["path_plot_network"]=str(folder_out+"plot_network.png")
metadata["saved_model"] = False
metadata["path_folder_save"]=str(folder_out)
metadata["suffix_write"]=str("train-"      + str(metadata["train_set_size"])              + "_" +
                             "test-"       + str(metadata["test_set_size"])               + "_" +
                             "layers-"     + str(metadata["n_hidden_layer"])              + "_" +
                             "n_nodes-"    + str(metadata["n_nodes_per_layer"])           + "_" + 
                             "nmaxSOAP-"   + str(metadata["nmax_SOAP"])                   + "_" + 
                             "lmaxSOAP-"   + str(metadata["lmax_SOAP"])                   + "_" +
                             "sigmaSOAP-"  + str(metadata["sigma_SOAP"])                  + "_" + 
                             "cutoffSOAP-" + str(metadata["cutoff_SOAP"])                 + "_" +
                             "drop_out0-"  + str(metadata["dropout_coef"][0,0])           + "_" + 
                             "drop_outN-"  + str(metadata["dropout_coef"][1,0])           ) 


model, metadata_stat, predictions_train, predictions_test = behler.buildTrainPredictWrite( input_train,input_test,output_train,output_test, metadata["species"], metadata["n_species"], metadata["n_features"], metadata["start_species"], metadata["nb_element_species"], metadata["n_nodes_structure"], metadata["dropout_coef"], metadata["kernel_constraint"],
                           activation_fct = metadata["activation_fct"],
                           loss_fct=metadata["loss_fct"],
                           optimizer=metadata["optimizer"],
                           path_folder_save=metadata["path_folder_save"], 
                           early_stop_metric=metadata["early_stop_metric"],
                           plot_network=metadata["plot_network"],
                           path_plot_network=metadata["path_plot_network"],
                           suffix_write=metadata["suffix_write"])

# Write the comparative between predictions and outputs
file_comp_train = str( metadata["path_folder_save"] + "ComparativeErrorsTrain_" +metadata["suffix_write"] )
file_comp_test  = str( metadata["path_folder_save"] + "ComparativeErrorsTest_"  +metadata["suffix_write"] )
metadata, output_train = mtd.deScaleEnergy( output_train_scale, metadata )
metadata, output_test  = mtd.deScaleEnergy( output_test_scale,  metadata )
metadata, predictions_train = mtd.deScaleEnergy( predictions_train, metadata )
metadata, predictions_test  = mtd.deScaleEnergy( predictions_test,  metadata )
behler.writeComparativePrediction(file_comp_train, output_train, predictions_train )
behler.writeComparativePrediction(file_comp_test,  output_test, predictions_test   )



energies = mtd.deScaleEnergy( behler.getAtomicEnergy( metadata["species"][0], 
                                   metadata["n_nodes_structure"][0,:],
                                   metadata["dropout_coef"][0,:],
                                   metadata["activation_fct"],
                                   metadata["loss_fct"], 
                                   metadata["optimizer"], 
                                   metadata["kernel_constraints"], 
                                   metadata["early_stop_metrics"],
                                   input_test_scale, 
                                   model ))

