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
import cpmd
import filexyz

from sklearn.preprocessing import StandardScaler  
from sklearn.decomposition import PCA


#==============================================================================
def handleOptionArg( input_label, default_value, metadata ):
    if not input_label in metadata :
        metadata[input_label] = default_value
    return metadata
def handleAllOptionaArg( input_label_list, default_values_list, metadata ):
    for i in range( len(input_label_list) ):
        metadata = handleOptionArg( input_label_list[i], default_values_list[i], metadata )
    return metadata
#==============================================================================

def getIndexExcluding( index_excluding, size_total ):
    return np.array(list(filter(lambda x : x not in index_excluding, np.arange(size_total))))

#============================================================================== 
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
    return metadata, structures_train, energies_train 
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
def choseTestDataByIndex(metadata,structures,energies,chosen_index): 
    metadata['test_index'] = chosen_index
    if metadata['test_set_size'] == 0 :    
        metadata['test_set_size']=len(chosen_index)
    structures_test = np.empty(metadata['test_set_size'],dtype=ase.atoms.Atoms)
    energies_test = np.empty(metadata['test_set_size'],dtype=float)
    if type(structures[0]) != ase.atoms.Atoms :
        for i in range(metadata['test_set_size']) :
            energies_test[i] = energies[chosen_index[i]]
            structures_test[i] = ase.atoms.Atoms(numbers=metadata['n_atoms'], positions=structures[chosen_index[i],:] )   
    else:
        for i in range(metadata['test_set_size']):
            energies_test[i] = energies[chosen_index[i]]          
            structures_test[i] = structures[chosen_index[i]]
    return metadata, structures_test, energies_test
#------------------------------------------------------------------------------
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
#==============================================================================

# Principal Component Analysis stuff
#==============================================================================
def makePCA( n_components ):
    return PCA(n_components=n_components)
#------------------------------------------------------------------------------
default_path_pcavar = "./"
def pcaVariance( input_, path_pcavar=default_path_pcavar ):
    nb_point=len(input_[:])*len(input_[0])
    pca=PCA().fit( np.array(input_[:]).reshape(nb_point,np.shape(input_[0])[1]) )
    var=np.cumsum(pca.explained_variance_ratio_)
    file_out=open(path_pcavar,"w")
    for i in range( np.shape(var)[0] ):
        file_out.write( str(i)+" "+str(var[i])+"\n" )
    file_out.close()
    return var
#------------------------------------------------------------------------------
def pcaSelectBestParams(descriptors,meta):
    pca=[]
    for specie in range(len(metadata["species"])):
        pca.append(PCA(n_components=metadata["pca_n"]).fit(descriptors[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].reshape(descriptors[:,metadata["start_species"][specie]:metadata["start_species"][specie]+metadata["nb_element_species"][specie],:].shape[0]*metadata["nb_element_species"][specie],metadata["n_features"])))
        if metadata['verbose']: 
            print("Precision of new features ",metadata["species"][specie]," :",np.cumsum(pca[specie].explained_variance_ratio_)[-1],"\n" )        
    return pca
#==============================================================================

# Scalers - to be generalized
#==============================================================================
def createScaler( input_, metadata ):
    # Create a scaler for the input using Scitkit Learn Standard Scaler
    scalers = []
    for specie in range(metadata["n_species"]):
        scalers.append(StandardScaler())
        nb_atoms_total=metadata["nb_element_species"][specie]*metadata["train_set_size"]
        start_specie = metadata["start_species"][specie]
        end_specie   = start_specie + metadata["nb_element_species"][specie]
        scalers[specie].fit(np.array(input_[start_specie:end_specie]).reshape(nb_atoms_total,metadata['n_features']))    
    return scalers
#------------------------------------------------------------------------------
def applyScale( scalers, input_, metadata ):
    # Applies a created scaler
    input_array=np.array(input_[:])
    for specie in range(metadata["n_species"]):
        start_specie = metadata["start_species"][specie]
        end_specie   = start_specie + metadata["nb_element_species"][specie]
        nb_points = np.shape(input_array[:,:,:])[1]
        nb_atoms_total = metadata["nb_element_species"][specie]*nb_points
        input_array[start_specie:end_specie,:,:] = scalers[specie].transform( input_array[start_specie:end_specie,:,:].reshape(nb_atoms_total,metadata["n_features"]) ).reshape( metadata["nb_element_species"][specie], nb_points, metadata["n_features"] )
    for i in range(np.shape(input_array)[0]):
        input_[i] = input_array[i,:,:]
    return input_
#------------------------------------------------------------------------------
#  Basic hand_made scalers
def scaleData( data ):
    min_data = data.min()
    range_data = data.max() - min_data
    return ( data - min_data)/range_data, min_data, range_data
#------------------------------------------------------------------------------
def deScaleData( data_scaled, min_data, range_data ):
    return data_scaled*range_data+min_data
#==============================================================================


# Handling Atoms stuff - may be moved at a later point
#==============================================================================
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
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
def getNbAtomsPerSpecies( atoms, metadata):
    if type(atoms) == list:
        atoms=atoms[0]
    metadata["nb_element_species"]=np.zeros(metadata["n_species"],dtype=int)
    for specie in range( metadata["n_species"] ):
        for atom in range( len(atoms) ):
            if metadata["species"][specie] == pT.z2Names(atoms.numbers[atom]) :
                metadata["nb_element_species"][specie] += 1
    return metadata
#------------------------------------------------------------------------------
def sortAtomsUniq( atoms ):
    # Sorting by decreasing Z 
    # BEHOLD THE BUBLE SORT OF DEATH - Needs to be changed asap by something 
    # less ineficient - probably are python functions that can do stuff
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
#------------------------------------------------------------------------------
def buildInput(data):
    input_=[]
    data_=np.stack(data["descriptor"].str[:].values)
    nb_atoms=np.shape(data_)[1]
    for atom in range(nb_atoms):
        input_.append(data_[:,atom,:])
    return input_
#------------------------------------------------------------------------------
def sortAtomsSpecie( atoms ) :
    if type(atoms) == list:
        for index in range(len(atoms)):
            atoms[index] = sortAtomsUniq(atoms[index])
        return atoms
    else:
        return sortAtomsUniq(atoms)
#------------------------------------------------------------------------------
def getStartSpecies( atoms, metadata ):
    if not metadata["species_sorted"] :
        atoms=sortAtomsSpecie(atoms)
        metadata["species_sorted"] = True
    if type(atoms) == list: 
        atoms=atoms[0]
    metadata["start_species"] = np.zeros(metadata["n_species"],dtype=int)
    for specie in range( metadata["n_species"] ):
        for atom in range( len(atoms) ):
            if metadata["species"][specie] == pT.z2Names(atoms.numbers[atom]) :
                metadata["start_species"][specie] = atom
                break
    return metadata
#==============================================================================

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
default_replace = False   # Whether or not we can select several times the same data point in the training set

# Default Descriptors parameters
default_descriptor   = 'None'

# Data treatment
default_PCA = False
default_pca_N = 0 
default_scale = False

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
                  replace = default_replace,
                  descriptor=default_descriptor, 
 #                 neigh_lim=default_neigh_lim,
                  pca=default_PCA,
                  pca_n=default_pca_N,
                  scale=default_scale, 
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
            'replace': replace,

             #Descriptor
            'descriptor_type': descriptor,    #Choice of descriptor
            'n_features': 1,
            
            # PCA or not PCA
            'pca_check': pca,             #False if no PCA, a number if PCA to select N first axis
            'pca_n': pca_n,
            'scaler': scale,             #False if no scaling, None if scaling
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