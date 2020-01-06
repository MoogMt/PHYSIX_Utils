#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 12:15:45 2019

@author:  julienh with modification from moogmt

Contains functions for handling data
"""

import numpy as np
import pandas as pd
import ase

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

def scaleDescriptors(data,descriptors,metadata):
    scaler = []
    for specie in range(len(metadata["species"])):
        scaler.append(StandardScaler())        
        scaler[specie].fit(descriptors[:,0:2,:].reshape(int(metadata['train_set_size']*metadata['n_specie'][specie]),metadata['n_features']))    
        descriptors[:,0:2,:] = scaler[specie].transform(descriptors[:,0:2,:].reshape(descriptors[:,0:2,:].shape[0]*2,metadata["n_features"])).reshape(descriptors.shape[0],2,metadata["n_features"])
    return descriptors, scaler

def pcaSelectBestParams(descriptors,metadata):
    pca=[]
    for specie in range(len(metadata["species"])):
        pca.append(PCA(n_components=metadata["pca_n"]).fit(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,metadata["n_features"])))
        if metadata['verbose']: 
            print("Precision of new features ",metadata["species"][specie]," :",np.cumsum(pca[specie].explained_variance_ratio_)[-1],"\n" )        
    return pca
    