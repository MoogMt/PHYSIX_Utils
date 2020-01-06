#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 10:11:28 2019

@author: moogmt
"""

import os
import ase
import numpy as np
import periodicTable as pT

def getNbAtoms( atoms ):
    if type(atoms) == ase.atoms.Atoms :
        return len(atoms)
    elif type(atoms) == np.ndarray :
        return len(atoms)
    elif type(atoms) == list :
        if type(atoms[0]) == ase.atoms.Atoms :
            return len(atoms[0])
        elif type(atoms[0]) == np.ndarray :
            return len(atoms)

def getSpecies( atoms ):
    if type(atoms) == list:
        atoms=atoms[0]
    types_=[]
    types_=np.append(types_,pT.z2Names(atoms.numbers[0]))
    for atom in range(getNbAtoms(atoms)):
        check = True
        for type_ in range(len(types_)):
            if pT.z2Names(atoms.numbers[atom]) == types_[type_]:
                check = False
        if check :
            types_=np.append(types_,pT.z2Names(atoms.numbers[atom]))
    return types_

def getNbAtomsPerSpecies( atoms, metadata):
    if type(atoms) == list:
        atoms=atoms[0]
    metadata["nb_per_species"]=np.zeros(( len(metadata["n_species"]) ))
    for specie in range( len(metadata["n_species"]) ):
        for atom in range( len(atoms) ):
            if metadata["species"][specie] == pT.names2Z(atoms.numbers[atom]) :
                metadata["nb_per_species"][specie] += 1
    return metadata

def getStartSpecies( atoms, metadata ):
    
    return metadata

# Default physics params
default_n_atoms        = 1
default_n_species      = 1
default_species        = []
default_masses         = np.zeros(default_n_atoms)
default_pbc            = None
default_periodic       = False
default_total_size_set = 0

# Default Descriptors parameters
default_descriptor   = 'None'

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

# Verbose
default_verbose = False
default_n_jobs=1


def buildMetaData( traj_file, energy_file, output_folder, 
                  n_atoms=default_n_atoms,
                  n_species=default_n_species,
                  species=default_species,
                  masses=default_masses,
                  pbc=default_pbc, 
                  periodic=default_periodic,
                  total_size_set = default_total_size_set,
                  descriptor=default_descriptor, 
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
                  suffix=default_suffix,
                  verbose=default_verbose,
                  n_jobs=default_n_jobs,
                  ):
    
    # Building Metadata Dic
    metadata={
            # Physics
            'n_atoms':n_atoms,
            'n_species':n_species,
            'species':species,
            'masses': masses,   
            'pbc': pbc,
            'periodic': periodic,
            'total_size_set': total_size_set,

             #Descriptor
            'descriptor_type': descriptor,    #Choice of descriptor
            'n_features': 1,
            
            # PCA or not PCA
            'pca_check': pca,             #False if no PCA, a number if PCA to select N first axis
            'pca_n': pca_n,
            'scaler': scale,             #False if no scaling, None if scaling
                        
            #Neural Net
            'train_fraction': train_fraction,
            'train_set_size': train_set_size,
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

            # Verbose
            'verbose': verbose,
            'n_jobs': n_jobs,

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