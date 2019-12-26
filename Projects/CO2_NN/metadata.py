#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 10:11:28 2019

@author: moogmt
"""

import sys
import os
import numpy as np

# Default physics params
default_n_atoms = 1
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
default_test_fraction=0.1 # Fraction of total dataset
default_replace = False
default_activation_fct = 'tanh'
default_loss_fct = 'mean_squared_error'
default_opt = 'adam'
default_n_epoch = 1000 # Number of epoch for optimization?
default_patience = 100 # ??
default_n_nodes_per_layer= 80
default_n_hidden_layer=2
default_n_nodes_structure=np.ones((default_n_atoms,default_n_hidden_layer))*default_n_nodes_per_layer

# Default I/O settingss
default_path_to_import_model = ""
default_prefix = ""


def buildMetaData( traj_file, energy_file, output_folder, 
                  pbc=default_pbc, 
                  masses=default_masses, 
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
                  path_to_import_model=default_path_to_import_model,
                  prefix=default_prefix
                  ):
    
    # Building Metadata Dic
    metadata={
            
            # Physics
            'pbc': pbc,
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

            # I/O        
            'traj_file': traj_file,
            'energy_file': energy_file,
            "path_2_model": path_to_import_model,    # To import from precedent NN, give datetime, or None
            'prefix_files': prefix,
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