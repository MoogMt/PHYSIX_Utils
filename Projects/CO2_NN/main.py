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

# Equilibrate sampling
#=============================================================================#
# Creating training set
n_jobs = 8 # Number of parallel cores to use (CPU)
fraction_choice=0.4
nb_bins = 10
comp_time_lim = 10
index_train = mtd.samplingExtremitiesFraction( fraction_choice, energies, comput_time, nb_bins, comp_time_lim  )
train_set_size = len(index_train)
input_train  = mtd.extractTrajectory( traj, index_train )
output_train = energies[ index_train ]
# Creating testing set
replace_elements=True
test_set_size = 1000
input_test, output_test = mtd.choseDataRandomExclusion( traj, energies, test_set_size, index_train, replace_elements )
# Clearing traj from memory 
traj=[]
#==============================================================================

# CREATING DESCRIPTORS
#==============================================================================
# Build descriptors from positions (train set only)
sigma_  = 0.9  # 3*sigma ~ 2.7A relatively large spread
cutoff_ = 6.0 # cut_off SOAP, 
nmax_   = 3
lmax_   = 3
# Train set
#-----------------------------------------------------------------------------
input_train  = desc.createDescriptorsAllSOAP( input_train, species, sigma_, cutoff_, nmax_, lmax_, periodic )
# Test set
#------------------------------------------------------------------------------
input_test   = desc.createDescriptorsAllSOAP( input_test, species, sigma_, cutoff_, nmax_, lmax_, periodic )
#------------------------------------------------------------------------------
n_features = np.shape(input_train)[2]
# Scaling 
#------------------------------------------------------------------------------
# Scaling Energy
output_train, min_output_train, range_output_train = mtd.scaleData( output_train )
output_test,  min_output_test,  range_output_test  = mtd.scaleData( output_test )
# Scaling Input
scalers = mtd.createScaler( input_train, species, start_species, nb_element_species, train_set_size, n_features ) # Create scaler on training set
input_train = mtd.applyScale( scalers, input_train, species, start_species, nb_element_species, n_features )
input_test  = mtd.applyScale( scalers, input_test,  species, start_species, nb_element_species, n_features )
#=============================================================================#

# BUILDING NETWORK
#=============================================================================#
# Parameters of the Neural net

#from keras.constraints import max_norm

# Iteration parameters
loss_fct = 'mean_squared_error' # Loss function in the NN
optimizer = 'Adam'                    # Choice of optimizers for training of the NN weights 
learning_rate = 0.0001
n_epochs = 2000                  # Number of epoch for optimization?
patience = 20                  # Patience for convergence
restore_weights = True
batch_size = 16
verbose_train = 1
early_stop_metric=['mse']

# Subnetorks structure
activation_fct = 'relu'  # Activation function in the dense hidden layers
n_nodes_per_layer = 40           # Number of nodes per hidden layer
n_hidden_layer = 5               # Number of hidden layers
n_nodes_structure=np.ones((n_species,n_hidden_layer),dtype=int)*n_nodes_per_layer # Structure of the NNs (overrides the two precedent ones)
kernel_constraint = None
bias_constraint = None

# Dropout rates
dropout_rate=np.zeros(( n_species, n_hidden_layer+1)) # Dropout for faster convergence (can be desactivated) 
dropout_rate[0,:]=0.2    # Drop out rate between initial descriptor and specie sub_network
dropout_rate[1:,:]=0.5   # Drop out rate inside the nodes of the specie sub_network
    
# Plot network
plot_network=True
path_plot_network=str(folder_out+"plot_network.png")
saved_model = False
path_folder_save=str(folder_out)
suffix_write="Otters_Test"

import behler

model, metadata_stat, predictions_train, predictions_test = behler.buildTrainPredictWrite( input_train, 
                                                                                          input_test,
                                                                                          output_train,
                                                                                          output_test, 
                                                                                          species, 
                                                                                          n_species, 
                                                                                          n_features, 
                                                                                          start_species, 
                                                                                          nb_element_species, 
                                                                                          n_nodes_structure, 
                                                                                          dropout_rate, 
                                                                                          n_epochs,
                                                                                          batch_size,
                                                                                          kernel_constraint=kernel_constraint,
                                                                                          bias_constraint=bias_constraint,
                                                                                          patience = patience,
                                                                                          activation_fct = activation_fct,
                                                                                          loss_fct=loss_fct,
                                                                                          optimizer=optimizer,
                                                                                          learning_rate=learning_rate,
                                                                                          path_folder_save=path_folder_save, 
                                                                                          early_stop_metric=early_stop_metric,
                                                                                          plot_network=plot_network,
                                                                                          path_plot_network=path_plot_network,
                                                                                          suffix_write=suffix_write)

# Descaling energies
output_train = mtd.deScaleData( output_train, min_output_train, range_output_train )
output_test  = mtd.deScaleData( output_test,  min_output_test,  range_output_test )
predictions_train = mtd.deScaleData( predictions_train, min_output_train, range_output_train )
predictions_test  = mtd.deScaleData( predictions_test,  min_output_test, range_output_test )

# Write the comparative between predictions and outputs
file_comp_train = str( path_folder_save + "ComparativeErrorsTrain_" + suffix_write )
file_comp_test  = str( path_folder_save + "ComparativeErrorsTest_"  + suffix_write )
behler.writeComparativePrediction( file_comp_train, output_train, predictions_train )
behler.writeComparativePrediction( file_comp_test,  output_test, predictions_test   )
#
#import matplotlib.pyplot as plt
#
#n_figure=1
#
#plt.figure(n_figure)
#plt.xlabel("E_{output} (Ry)")
#plt.ylabel("E_{prediction} (Ry)")
#plt.plot(output_train-output_train.min(), predictions_train-output_train.min(),"r.")
#plt.plot(output_test-output_test.min(), predictions_test-output_test.min(), "b.")
#plt.legend(["Train","Test"])
#plt.show()
#n_figure +=1
#
#plt.figure(n_figure)
#plt.xlabel("Structure Index (#)")
#plt.ylabel("E_{output}^{train}-E_{prediction}^{train} (Ry/CO_{2})")
#plt.plot( (output_train-predictions_train)/96*13.6,"r-")
#plt.show()
#n_figure +=1
#
#plt.figure(n_figure)
#plt.xlabel("Structure Index (#)")
#plt.ylabel("E_{output}^{train}-E_{prediction}^{train} (Ry/CO_{2})")
#plt.plot( (output_test-predictions_test)/96*13.6,"r-")
#plt.show()
#n_figure +=1
#
#plt.figure(n_figure)
#plt.xlabel("Structure Index (#)")
#plt.ylabel("Energy (Ry)")
#plt.plot(output_test,"r-")
#plt.plot(predictions_test, "b-")
#plt.legend(["Output (test)","Prediction (test)"])
#plt.show()
#n_figure +=1
#
#plt.figure(n_figure)
#plt.xlabel("Structure Index (#)")
#plt.ylabel("Energy (Ry)")
#plt.plot(output_train,"r-")
#plt.plot(predictions_train, "b-")
#plt.legend(["Output (train)","Prediction (train)"])
#plt.show()
#n_figure +=1
#
#energies_train = []
#energies_test  = []
#for specie in range( n_species ):
#    energies_train.append(mtd.deScaleData( behler.getAtomicEnergy( species[specie], 
#                                                                     start_species[specie],
#                                                                     nb_element_species[specie],
#                                                                     n_nodes_structure[specie,:],
#                                                                     dropout_rate[specie,:],
#                                                                     input_train_scale, 
#                                                                     model,
#                                                                     activation_fct=activation_fct,
#                                                                     loss_fct=loss_fct, 
#                                                                     optimizer=optimizer, 
#                                                                     kernel_constraint=kernel_constraint, 
#                                                                     early_stop_metric=early_stop_metric),
#                                                                     min_output_train, range_output_train ) )
#    energies_test.append(mtd.deScaleData( behler.getAtomicEnergy( species[specie], 
#                                                                    start_species[specie],
#                                                                    nb_element_species[specie],
#                                                                    n_nodes_structure[specie,:],
#                                                                    dropout_rate[specie,:],
#                                                                    input_test_scale, 
#                                                                    model,
#                                                                    activation_fct=activation_fct,
#                                                                    loss_fct=loss_fct, 
#                                                                    optimizer=optimizer, 
#                                                                    kernel_constraint=kernel_constraint, 
#                                                                    early_stop_metric=early_stop_metric),
#                                                                    min_output_test, range_output_test ) )
#
#
#nb_bins=100
#plt.figure(1)
#for specie in range( n_species ):
#    plt.xlabel("Energy (Ry)")
#    plt.hist(energies_train[specie], bins=nb_bins)
#plt.show()
#n_figure+=1 
#
#plt.figure(2)
#for specie in range( n_species ):
#    plt.ylabel("Energy (Ry)")
#    plt.hist(energies_test[specie], bins=nb_bins)
#plt.show()
#n_figure+=1
