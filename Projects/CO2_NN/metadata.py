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
default_NN=np.zeros((1,1))

# Neural Nets default parameters
default_n_nodes=8
default_test_fraction=0.1 # Fraction of total dataset
default_activation_fct = 'tanh'
default_loss_fct = 'mean_squared_error'
default_opt = 'adam'
default_n_epoch = 1000 # Number of epoch for optimization?
default_patience = 100 # ??

# Default I/O settings
default_replace = False
default_import = None
default_prefix = ""


def buildMetaData( traj_file, energy_file, output_folder, pbc=default_pbc ):
    
    # Building Metadata Dic
    metadata={
            
            # Physics
            'pbc': pbc,
            'masses': np.zeros(default_n_atoms),   

             #Descriptor
            'descriptor_type': default_descriptor,    #Choice of descriptor
            #SOAP
            'scaler': default_scale,             #False if no scaling, None if scaling
            'sigma_SOAP':default_sigma_SOAP,
            'rcut': default_cutoff_SOAP, 
            'nmax': default_nmax_SOAP,
            'lmax': default_lmax_SOAP,
            # Local PIV
            'nearest_neigh':default_NN,
            # PCA or not PCA
            'N_PCA': default_PCA,             #False if no PCA, a number if PCA to select N first axis
                        
            #Neural Net
            'test_size':default_test_fraction,
            'activation_function': default_activation_fct,
            'loss_function': default_loss_fct,
            'optimizer': default_opt, 
            'epochs': default_n_epoch,
            'patience': default_patience,
            'N_nodes': np.ones((default_n_atoms,2))*default_n_nodes,
            'N_feature': 1,

            # I/O            
            'replace': default_replace,       # Choose whether you can pick twice the same data point in the train set at random
            "import_from": default_import,    # To import from precedent NN, give datetime, or None
            "output_folder": output_folder,
            'traj_file': traj_file,
            'energy_file': energy_file,
            'prefix_files': default_prefix,

            }
    return metadata

def checkMetaDataIO( metadata, verbose ):
    if not os.path.isfile( metadata['traj_file'] ): 
        if verbose:
            print("Necessary input trajectory file does not exists! Exiting...")
        return False
    elif not os.path.isfile( metadata['energy_file'] ):
        if verbose:
            print("Necessary input Energy file does not exists! Exiting")
        return False
    return True