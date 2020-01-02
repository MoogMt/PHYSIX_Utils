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
    if metadata['size_train_set'] == 0 :    
        metadata['size_train_set']=len(chosen_index)
    structures_train = np.empty(metadata['size_train_set'],dtype=ase.atoms.Atoms)
    energies_train = np.empty(metadata['size_train_set'],dtype=float)
    if type(structures[0]) != ase.atoms.Atoms :
        for i in range(metadata['size_train_set']) :
            energies_train[i] = energies[chosen_index[i]]
            structures_train[i] = ase.atoms.Atoms(numbers=metadata['n_atoms'], positions=structures[chosen_index[i],:] )   
    else:
        for i in range(metadata['size_train_set']):
            energies_train[i] = energies[chosen_index[i]]          
            structures_train[i] = structures[chosen_index[i]]
    return metadata, pd.DataFrame({"energy":energies_train,'structures':structures_train})

# 
def choseTrainDataRandom(metadata,structures,energies):    
    chosen_index = np.random.choice(metadata['total_size_set'],size=int(metadata['train_fraction']*metadata['total_size_set']),replace=metadata['replace'])
    return choseTrainDataByIndex(metadata,structures,energies,chosen_index)

def pcaSelectBestParams(descriptors,nb_features,N_PCA):
    pca_O = PCA(n_components=N_PCA).fit(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features))
    pca_H = PCA(n_components=N_PCA).fit(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features))
    print("Precision of new features [O,H] = ",[np.cumsum(pca_O.explained_variance_ratio_)[-1],np.cumsum(pca_H.explained_variance_ratio_)[-1]])        
    return [pca_O, pca_H]
    
def scale_descriptors(data,descriptors):
    is_train = data['is_train']
    nb_features = descriptors.shape[2]
    scaler = []
    scaler.append(StandardScaler())
    scaler[0].fit(descriptors[is_train][:,0:2,:].reshape(int(is_train.sum()*2),nb_features))    
    descriptors[:,0:2,:] = scaler[0].transform(descriptors[:,0:2,:].reshape(descriptors[:,0:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,nb_features)
    scaler.append(StandardScaler())
    scaler[1].fit(descriptors[is_train][:,2:,:].reshape(int(is_train.sum()*5),nb_features))    
    descriptors[:,2:,:] = scaler[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,nb_features)
    return descriptors, scaler