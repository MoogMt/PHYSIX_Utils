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
data_out     = cpmd.readEnergiesFile( file_energies )
energies     = cpmd.extractPotentialEnergy(     data_out )
comput_time  = cpmd.extractSCFcomputationTime(  data_out )
# Homogenizing with stride
stride_energies=5 # Improve by reading the stride in the input file and making adjustements.
energies    = energies    [ 0: len(energies):stride_energies ]
comput_time = comput_time [ 0: len(energies):stride_energies ]
# add a check to verify congruence of sizes...
# Getting species present in the simulation
metadata['n_atoms'] = len(traj[0])
metadata['species'] = mtd.getSpecies(traj[0])
metadata['n_species'] = len(metadata['species'])
#traj=mtd.sortAtomsSpecie(traj) # Use if you need to sort the atoms per type, current implementation is *very* slow
metadata["species_sorted"]=True
metadata['total_size_set'] = len(energies)
metadata=mtd.getStartSpecies(traj,metadata)
metadata=mtd.getNbAtomsPerSpecies(traj,metadata)
#=============================================================================#

# CREATING DESCRIPTORS
#=============================================================================#
# Creating training set
n_jobs = 8 # Number of parallel cores to use (CPU)
fraction_choice=0.4
nb_bins = 10
comp_time_lim = 10
index_train = mtd.samplingExtremitiesFraction( fraction_choice, energies, comput_time, nb_bins, comp_time_lim  )
train_set_size = len(index_train)
input_train_raw  = mtd.extractTrajectory( traj, index_train )
output_train_raw = energies[ index_train ]
# Creating testing set
test_set_size = 200
metadata, input_test_raw, output_test_raw = mtd.choseTestDataRandomExclusion(metadata,traj,energies)
# Build descriptors from positions (train set only)
sigma_  = 0.9  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 4.0 # cut_off SOAP, 
nmax_   = 3
lmax_   = 2
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
output_train_scale, min_output_train, range_output_train = mtd.scaleData( output_train_raw, metadata )
output_test_scale,  min_output_test,  range_output_test  = mtd.scaleData( output_test_raw,  metadata )
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
metadata["patience"] = 20                  # Patience for convergence
metadata["restore_weights"] = True
metadata["batch_size"] = 300
metadata["verbose_train"] = 1
metadata["early_stop_metric"]=['mse']

# Subnetorks structure
metadata["activation_fct"] = 'relu'  # Activation function in the dense hidden layers
metadata["n_nodes_per_layer"] = metadata["n_features"]           # Number of nodes per hidden layer
metadata["n_hidden_layer"] = 3               # Number of hidden layers
metadata["n_nodes_structure"]=np.ones((metadata["n_species"],metadata["n_hidden_layer"]),dtype=int)*metadata["n_nodes_per_layer"] # Structure of the NNs (overrides the two precedent ones)
metadata["kernel_constraint"] = None

# Dropout rates
metadata["dropout_coef"]=np.zeros((metadata["n_species"],metadata["n_hidden_layer"]+1)) # Dropout for faster convergence (can be desactivated) 
metadata["dropout_coef"][0,:]=0.2    # Drop out rate between initial descriptor and specie sub_network
metadata["dropout_coef"][1:,:]=0.5   # Drop out rate inside the nodes of the specie sub_network
    
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


import behler

model, metadata_stat, predictions_train, predictions_test = behler.buildTrainPredictWrite( input_train_scale,input_test_scale,output_train_scale,output_test_scale, 
                           metadata["species"], 
                           metadata["n_species"], 
                           metadata["n_features"], 
                           metadata["start_species"], 
                           metadata["nb_element_species"], 
                           metadata["n_nodes_structure"], 
                           metadata["dropout_coef"], 
                           metadata["n_epochs"],
                           metadata["batch_size"],
                           metadata["kernel_constraint"],
                           activation_fct = metadata["activation_fct"],
                           loss_fct=metadata["loss_fct"],
                           optimizer=metadata["optimizer"],
                           path_folder_save=metadata["path_folder_save"], 
                           early_stop_metric=metadata["early_stop_metric"],
                           plot_network=metadata["plot_network"],
                           path_plot_network=metadata["path_plot_network"],
                           suffix_write=metadata["suffix_write"])

# Descaling energies
output_train = mtd.deScaleData( output_train_scale, min_output_train, range_output_train )
output_test  = mtd.deScaleData( output_test_scale,  min_output_test,  range_output_test )
predictions_train = mtd.deScaleData( predictions_train, metadata )
predictions_test  = mtd.deScaleData( predictions_test,  metadata )

# Write the comparative between predictions and outputs
file_comp_train = str( metadata["path_folder_save"] + "ComparativeErrorsTrain_" +metadata["suffix_write"] )
file_comp_test  = str( metadata["path_folder_save"] + "ComparativeErrorsTest_"  +metadata["suffix_write"] )
behler.writeComparativePrediction( file_comp_train, output_train, predictions_train )
behler.writeComparativePrediction( file_comp_test,  output_test, predictions_test   )

import matplotlib.pyplot as plt

n_figure=1

plt.figure(n_figure)
plt.xlabel("E_{output} (Ry)")
plt.ylabel("E_{prediction} (Ry)")
plt.plot(output_train-output_train.min(), predictions_train-output_train.min(),"r.")
plt.plot(output_test-output_test.min(), predictions_test-output_test.min(), "b.")
plt.legend(["Train","Test"])
plt.show()
n_figure +=1

plt.figure(n_figure)
plt.xlabel("Structure Index (#)")
plt.ylabel("E_{output}^{train}-E_{prediction}^{train} (Ry/CO_{2})")
plt.plot( (output_train-predictions_train)/96*13.6,"r-")
plt.show()
n_figure +=1

plt.figure(n_figure)
plt.xlabel("Structure Index (#)")
plt.ylabel("E_{output}^{train}-E_{prediction}^{train} (Ry/CO_{2})")
plt.plot( (output_test-predictions_test)/96*13.6,"r-")
plt.show()
n_figure +=1

plt.figure(n_figure)
plt.xlabel("Structure Index (#)")
plt.ylabel("Energy (Ry)")
plt.plot(output_test,"r-")
plt.plot(predictions_test, "b-")
plt.legend(["Output (test)","Prediction (test)"])
plt.show()
n_figure +=1

plt.figure(n_figure)
plt.xlabel("Structure Index (#)")
plt.ylabel("Energy (Ry)")
plt.plot(output_train,"r-")
plt.plot(predictions_train, "b-")
plt.legend(["Output (train)","Prediction (train)"])
plt.show()
n_figure +=1

energies_train = []
energies_test  = []
for specie in range( metadata["n_species"] ):
    energies_train.append(mtd.deScaleEnergy( behler.getAtomicEnergy( metadata["species"][specie], 
                                                                   metadata["start_species"][specie],
                                                                   metadata["nb_element_species"][specie],
                                                                   metadata["n_nodes_structure"][specie,:],
                                                                   metadata["dropout_coef"][specie,:],
                                                                   input_train_scale, 
                                                                   model,
                                                                   activation_fct=metadata["activation_fct"],
                                                                   loss_fct=metadata["loss_fct"], 
                                                                   optimizer=metadata["optimizer"], 
                                                                   kernel_constraint=metadata["kernel_constraint"], 
                                                                   early_stop_metric=metadata["early_stop_metric"]),
                                                                    metadata )[1])
    energies_test.append(mtd.deScaleEnergy( behler.getAtomicEnergy( metadata["species"][specie], 
                                                                   metadata["start_species"][specie],
                                                                   metadata["nb_element_species"][specie],
                                                                   metadata["n_nodes_structure"][specie,:],
                                                                   metadata["dropout_coef"][specie,:],
                                                                   input_test_scale, 
                                                                   model,
                                                                   activation_fct=metadata["activation_fct"],
                                                                   loss_fct=metadata["loss_fct"], 
                                                                   optimizer=metadata["optimizer"], 
                                                                   kernel_constraint=metadata["kernel_constraint"], 
                                                                   early_stop_metric=metadata["early_stop_metric"]),
                                                                    metadata )[1])


nb_bins=100
plt.figure(1)
for specie in range(metadata["n_species"]):
    plt.hist(energies_train[specie], bins=nb_bins)
plt.show()
n_figure+=1 

plt.figure(2)
for specie in range(metadata["n_species"]):
    plt.hist(energies_test[specie], bins=nb_bins)
plt.show()
n_figure+=1