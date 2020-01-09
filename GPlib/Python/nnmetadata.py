#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 12:15:45 2019

@author:  julienh with modification from moogmt

Contains functions for handling atomic simulation data
"""

import os
import ase
import numpy as np
import periodicTable as pT
import pandas as pd

from sklearn.preprocessing import StandardScaler  
from sklearn.decomposition import PCA

# 
def choseTrainDataByIndex(metadata,structures,energies,chosen_index): 
    metadata['train_index'] = chosen_index
    if metadata['train_set_size'] == 0 :    
        metadata['train_set_size']=len(chosen_index)
    structures_train = np.empty(metadata['train_set_size'],dtype=ase.atoms.Atoms)
    energies_train = np.empty(metadata['train_set_size'],dtype=float)
    if type(structures[0]) != ase.atoms.Atoms :
        for i in range(metadata['train_set_size']) :
            energies_train[i] = energies[chosen_index[i]]
            structures_train[i] = ase.atoms.Atoms(numbers=metadata['n_atoms'], positions=structures[chosen_index[i],:] )   
    else:
        for i in range(metadata['train_set_size']):
            energies_train[i] = energies[chosen_index[i]]          
            structures_train[i] = structures[chosen_index[i]]
    return metadata, pd.DataFrame({"energy":energies_train,'structures':structures_train})

# 
def choseTrainDataRandom(metadata,structures,energies):    
    if metadata['train_set_size'] == 0 :
        if metadata['train_fraction'] < 1:
            print("Invalid train_faction in metadata","\n")
            return False, False
        if metadata['total_size_set'] < 1 :
            print("Invalid total_size_set in metadata","\n")
        metadata['train_set_size'] = int(metadata['train_fraction']*metadata['total_size_set'])
    chosen_index = np.random.choice(metadata['total_size_set'],size=metadata['train_set_size'],replace=metadata['replace'])
    return choseTrainDataByIndex(metadata,structures,energies,chosen_index)

#  THIS SHOULD BE IMPROVED **SOMEWHAT**
def choseTestDataByIndex(metadata,structures,energies,chosen_index): 
    metadata['test_index'] = chosen_index
    if metadata['test_set_size'] == 0 :    
        metadata['test_set_size']=len(chosen_index)
    structures_train = np.empty(metadata['test_set_size'],dtype=ase.atoms.Atoms)
    energies_train = np.empty(metadata['test_set_size'],dtype=float)
    if type(structures[0]) != ase.atoms.Atoms :
        for i in range(metadata['test_set_size']) :
            energies_train[i] = energies[chosen_index[i]]
            structures_train[i] = ase.atoms.Atoms(numbers=metadata['n_atoms'], positions=structures[chosen_index[i],:] )   
    else:
        for i in range(metadata['test_set_size']):
            energies_train[i] = energies[chosen_index[i]]          
            structures_train[i] = structures[chosen_index[i]]
    return metadata, pd.DataFrame({"energy":energies_train,'structures':structures_train})

def getIndexExcluding( index_excluding, size_total ):
    return np.array(list(filter(lambda x : x not in index_excluding, np.arange(size_total))))

def choseTestDataRandomExclusion(metadata,structures,energies):    
    if metadata['test_set_size'] == 0 :
        if metadata['test_fraction'] < 1:
            print("Invalid train_faction in metadata","\n")
            return False, False
        if metadata['total_size_set'] < 1 :
            print("Invalid total_size_set in metadata","\n")
        metadata['train_set_size'] = int(metadata['test_fraction']*metadata['total_size_set'])
    chosen_index=getIndexExcluding( metadata["train_index"], metadata["total_size_set"] )
    return choseTestDataByIndex(metadata,structures,energies,chosen_index)

def scaleData(data,metadata):
    scaler = []
    for specie in range(metadata["n_species"]):
        scaler.append(StandardScaler())        
        scaler[specie].fit(data[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].reshape(metadata['train_set_size']*metadata['n_species'][specie],metadata['n_features']))    
        data[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:] = scaler[specie].transform(data[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].reshape(data[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].shape[0]*metadata["nb_element_species"][specie],metadata["n_features"])).reshape(data.shape[0],metadata["nb_element_species"][specie],metadata["n_features"])
    return data, scaler

def pcaSelectBestParams(descriptors,metadata):
    pca=[]
    for specie in range(len(metadata["species"])):
        pca.append(PCA(n_components=metadata["pca_n"]).fit(descriptors[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].reshape(descriptors[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].shape[0]*metadata["nb_element_species"][specie],metadata["n_features"])))
        if metadata['verbose']: 
            print("Precision of new features ",metadata["species"][specie]," :",np.cumsum(pca[specie].explained_variance_ratio_)[-1],"\n" )        
    return pca

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
    metadata["nb_element_species"]=np.zeros(metadata["n_species"],dtype=int)
    for specie in range( metadata["n_species"] ):
        for atom in range( len(atoms) ):
            if metadata["species"][specie] == pT.names2Z(atoms.numbers[atom]) :
                metadata["nb_element_species"][specie] += 1
    return metadata

def sortAtomsUniq( atoms ):
    # Sorting by decreasing Z 
    # BEHOLD THE BUBLE SORT OF DEATH - Y a probablement plus efficace
    # Mais sachant qu'on fera *jamais* plus de tableaux de 10^4 elements,
    # et que c'est pas qqchose qu'on va faire régulièrement, ça va...
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if atoms.numbers[i] < atoms.numbers[j]:
                store_z = atoms.numbers[i]
                store_positions = atoms.positions[i,:]
                atoms.numbers[i] = atoms.numbers[j]
                atoms.positions[i,:] = atoms.positions[j,:]
                atoms.numbers[j] = store_z
                atoms.positions[j,:] = store_positions
    return atoms

def sortAtomsSpecie( atoms ) :
    if type(atoms) == list:
        for index in range(len(atoms)):
            atoms[index] = sortAtomsUniq(atoms[index])
        return atoms
    else:
        return sortAtomsUniq(atoms)

def getStartSpecies( atoms, metadata ):
    if not metadata["species_sorted"] :
        atoms=sortAtomsSpecie(atoms)
        metadata["species_sorted"] = True
    if type(atoms) == list: 
        atoms=atoms[0]
    metadata["start_species"] = np.zeros(metadata["n_species"],dtype=int)
    for specie in range( metadata["n_species"] ):
        for atom in range( len(atoms) ):
            if metadata["species"][specie] == pT.names2Z(atoms.numbers[atom]) :
                metadata["start_species"] = atom
                break
    return metadata

# Default physics params
default_n_atoms        = 1
default_n_species      = 1
default_species        = []
default_masses         = np.zeros(default_n_atoms)
default_pbc            = None
default_periodic       = False
default_total_size_set = 0
default_start_species  = np.zeros(default_n_atoms,dtype=int)
default_nb_element_species = np.zeros(default_n_atoms,dtype=int)
default_species_sorted = False

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
                  start_species=default_start_species,
                  nb_element_species=default_nb_element_species,
                  species_sorted=default_species_sorted,
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
            'start_species': start_species,
            'nb_element_species': nb_element_species,
            'species_sorted': species_sorted,

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