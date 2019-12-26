#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 10:11:28 2019

@author: moogmt
"""

import sys
import os
import numpy as np

# Default Descriptors parameters
default_descriptor = 'NN'
default_PCA = False
# - SOAP
default_scale = False
default_sigma_SOAP = 0.8
default_cutoff_SOAP = 3.5 # Angstroms
default_nmax_SOAP = 3
default_lmax_SOAP = 2
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

def buildMetaData( input_folder, output_folder, volume, temperature, n_atoms ):
    
    # Building Metadata Dic
    metadata={
            
            # Physics
            'temperature': temperature,
            'volume': volume,
            'masses': np.zeros(n_atoms),   

             #Descriptor
            'descriptor_type': default_descriptor,    #Choice of descriptor
            'N_PCA': default_PCA,             #False if no PCA, a number if PCA to select N first axis
            #SOAP
            'scaler': default_scale,             #False if no scaling, None if scaling
            'sigma_SOAP':default_sigma_SOAP,
            'rcut': default_cutoff_SOAP, 
            'nmax': default_nmax_SOAP,
            'lmax': default_lmax_SOAP,
                        
            #Neural Net
            'test_size':default_test_fraction,
            'activation_function': default_activation_fct,
            'loss_function': default_loss_fct,
            'optimizer': default_opt, 
            'epochs': default_n_epoch,
            'patience': default_patience,
            'N_nodes': np.ones((n_atoms,2))*default_n_nodes,

            # I/O            
            'replace': default_replace,
            "import_from": default_import,         #To import from precedent NN, give datetime, or None
            "output_folder": output_folder,
            'input_folder': input_folder,
            'prefix_files': default_prefix,

            }
    return metadata

def checkMetaDataIO( metadata, verbose ):
    if not os.path.isdir(metadata['input_folder']) :
        if verbose:
            print("Input folder does not exists!")
        return False
    elif not os.path.isdir(metadata['output_folder'] ) : 
        if verbose:
            print("Output folder does not exists, creating it.")
        os.mkdir(metadata['output_folder'])
    return True