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
default_total_size_set = 0

# Default Descriptors parameters
default_descriptor = 'SOAP'
# - SOAP
default_sigma_SOAP = 0.
default_cutoff_SOAP = 0. # Angstroms
default_nmax_SOAP = 0
default_lmax_SOAP = 0
# - NN
#default_neigh_lim=np.zeros((1,1))
# - ACSF
default_cutoff_ACSF = 0.                          # Cut-off for the ACSF function
default_n_acsf = 0                                # Number of ACSF functions per atoms
default_g2_params = np.zeros((default_n_acsf,2))  # Parameters for the g2 functions for the ACSF functions
default_g3_params = np.zeros((default_n_acsf,3))  # Parameters for the g3 functions for the ACSF functions

# Data treatment
default_PCA = False
default_pca_N = 0 
default_scale = False

# Neural Net default parameters
default_train_fraction = 0.1 # Fraction of total dataset, if chosing data randomly from total set
default_train_set_size = 0   # Size of the training set
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
                  total_size_set = default_total_size_set,
                  descriptor=default_descriptor, 
                  sigma_SOAP=default_sigma_SOAP, 
                  cutoff_SOAP=default_cutoff_SOAP, 
                  nmax_SOAP=default_nmax_SOAP, 
                  lmax_SOAP=default_lmax_SOAP,
                  cutoff_ACSF=default_cutoff_ACSF,
                  n_acsf=default_n_acsf,
                  g2_params=default_g2_params,
                  g3_params=default_g3_params,
 #                 neigh_lim=default_neigh_lim,
                  pca=default_PCA,
                  pca_n=default_pca_N,
                  scale=default_scale, 
                  train_fraction=default_train_fraction,
                  train_set_size=default_train_set_size,
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
            'n_atoms':n_atoms,
            'n_species':n_species,
            'masses': masses,   
            'pbc': pbc,
            'total_size_set': total_size_set,

             #Descriptor
            'descriptor_type': descriptor,    #Choice of descriptor
            'n_features': 1,
            # - SOAP
            'sigma_SOAP': sigma_SOAP,
            'rcut': cutoff_SOAP, 
            'nmax': nmax_SOAP,
            'lmax': lmax_SOAP,
            # - NN
 #           'nearest_neigh': neigh_lim,
            # - ACSF
            
            # PCA or not PCA
            'pca_check': pca,             #False if no PCA, a number if PCA to select N first axis
            'pca_n': pca_n,
            'scaler': scale,             #False if no scaling, None if scaling
                        
            #Neural Net
            'train_fraction': train_fraction,
            'size_train_set': train_set_size,
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