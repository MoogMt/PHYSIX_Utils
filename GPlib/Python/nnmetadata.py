#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 10:11:28 2019

@author: moogmt
"""

import os
import numpy as np

# Default physics params
default_n_atoms = 1
default_n_species = 1
default_masses=np.zeros(default_n_atoms)
default_pbc=None

# Default Descriptors parameters
default_descriptor = 'SOAP'
default_PCA = False
# - SOAP
default_scale = False
default_sigma_SOAP = 0.8
default_cutoff_SOAP = 3.5 # Angstroms
default_nmax_SOAP = 3
default_lmax_SOAP = 2
# - NN
default_neigh_lim=np.zeros((1,1))

# Neural Net default parameters
default_test_fraction=0.1 # Fraction of total dataset, if chosing data randomly from total set
default_replace = False   # Whether or not we can select several times the same data point in the training set
default_activation_fct = 'tanh'  # Activation function in the dense hidden layers
default_loss_fct = 'mean_squared_error' # Loss function in the NN
default_opt = 'adam'                    # Choice of optimizers for training of the NN weights 
default_n_epoch = 1000                  # Number of epoch for optimization?
default_patience = 100                  # Patience for convergence
default_n_nodes_per_layer= 80           # Number of nodes per hidden layer
default_n_hidden_layer=2                # Number of hidden layers
default_n_nodes_structure=np.ones((default_n_species,default_n_hidden_layer))*default_n_nodes_per_layer # Structure of the NNs (overrides the two precedent ones)
default_dropout_coef=np.zeros((default_n_hidden_layer+1,default_n_species)) # Dropout for faster convergence (can be desactivated) 

# Default I/O settingss
default_path_to_import_model = ""        # If taking over from previous model, path where to save the NN
default_prefix = ""                      # prefix for the output files
default_suffix = ""                      # suffix for the output files


def buildMetaData( traj_file, energy_file, output_folder, 
                  n_atoms=default_n_atoms,
                  n_species=default_n_species,
                  masses=default_masses, 
                  pbc=default_pbc, 
                  descriptor=default_descriptor, 
                  scale=default_scale, 
                  sigma_SOAP=default_sigma_SOAP, 
                  cutoff_SOAP=default_cutoff_SOAP, 
                  nmax_SOAP=default_nmax_SOAP, 
                  lmax_SOAP=default_lmax_SOAP,
                  neigh_lim=default_neigh_lim,
                  pCA=default_PCA,
                  test_fraction=default_test_fraction,
                  replace=default_replace,
                  activation_fct=default_activation_fct,
                  loss_fct=default_loss_fct,
                  opt=default_opt,
                  n_epoch=default_n_epoch,
                  patience=default_patience,
                  n_nodes_per_layer=default_n_nodes_per_layer,
                  n_hidden_layer=default_n_hidden_layer,
                  n_nodes_structure=default_n_nodes_structure,
                  dropout_coef=default_dropout_coef,
                  path_to_import_model=default_path_to_import_model,
                  prefix=default_prefix,
                  suffix=default_suffix
                  ):
    
    # Building Metadata Dic
    metadata={
            
            # Physics
            'pbc': pbc,
            'n_atoms':n_atoms,
            'n_species':n_species,
            'masses': masses,   

             #Descriptor
            'descriptor_type': descriptor,    #Choice of descriptor
            'n_features': 1,
            #SOAP
            'scaler': scale,             #False if no scaling, None if scaling
            'sigma_SOAP': sigma_SOAP,
            'rcut': cutoff_SOAP, 
            'nmax': nmax_SOAP,
            'lmax': lmax_SOAP,
            # Local PIV
            'nearest_neigh': neigh_lim,
            # PCA or not PCA
            'N_PCA': pCA,             #False if no PCA, a number if PCA to select N first axis
                        
            #Neural Net
            'test_size': test_fraction,
            'replace': replace,       # Choose whether you can pick twice the same data point in the train set at random
            'activation_function': activation_fct,
            'loss_function': loss_fct,
            'optimizer': opt, 
            'epochs': n_epoch,
            'patience': patience,
            'n_nodes_per_layer': n_nodes_per_layer,
            'n_hidden_layer': n_hidden_layer,
            'n_nodes_structure': n_nodes_structure,
            'dropout_coef': dropout_coef,

            # I/O        
            'traj_file': traj_file,
            'energy_file': energy_file,
            "path_2_model": path_to_import_model,    # To import from precedent NN, give datetime, or None
            'prefix_out': prefix,
            'suffix_out': suffix,
            "output_folder": output_folder,

            }
    return metadata

def checkMetaDataIO( metadata, verbose ):
    if not os.path.isfile( metadata['traj_file'] ): 
        if verbose:
            print("Necessary input trajectory file does not exists! Exiting...")
        return False
    elif not os.path.isfile( metadata['energy_file'] ):
        if verbose:
            print("Necessary input Energy file does not exists! Exiting...")
        return False
    if not os.path.isdir( metadata['output_folder']):
        if verbose:
            print("Output folder does not exists, creating it...")
        os.mkdir(metadata['output_folder'])
    return True