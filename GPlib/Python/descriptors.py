#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 18:04:53 2019

@author: moogmt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:49:53 2019

@author: julienh
"""
import numpy as np
import tqdm
import pandas as pd

from sklearn.preprocessing import StandardScaler  

from dscribe.descriptors import SOAP
from dscribe.descriptors import ACSF

from scipy.spatial.distance import cdist



def create_is_train(tot_time,test_size=0.15):
    is_train = np.ones(tot_time)
    is_train[:int(tot_time*test_size)] = 0
    is_train = np.random.choice(is_train,size=is_train.size,replace=False).astype(bool)
    return is_train

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

def create_data_SOAP(data, metadata):   
    # Prepping SOAP structure
    soap = SOAP( species=np.unique(metadata['n_atoms']), sigma=metadata['sigma_SOAP'], periodic=metadata['periodic'], rcut=metadata['cutoff_SOAP'], nmax=metadata['nmax'], lmax=metadata['lmax'],sparse=metadata['sparse_SOAP'])
    metadata['n_features'] = soap.get_number_of_features()
    
    # Creating descriptors
    descriptors = pd.np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_features']))
    for index in tqdm.tqdm(range(metadata['total_set_size'])):
        descriptors[index] = soap.create(data['structures'][index],positions=np.arange(metadata['n_atoms']))
    
    #selecting best params
    if metadata['n_PCA']:
            try:
                metadata['PCAs']
            except KeyError:
                PCAs = select_best_params(descriptors[data['is_train'].values],nb_features,N_PCA) 
                new_descriptors = pd.np.empty((data.index.max()+1,len(particles),N_PCA))
                new_descriptors[:,:2,:] = PCAs[0].transform(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,N_PCA)
                new_descriptors[:,2:,:] = PCAs[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,N_PCA)
                descriptors = new_descriptors
                metadata['old_N_feature'] = nb_features
                nb_features = N_PCA
                metadata['PCAs'] = PCAs
                

            else:
                PCAs = metadata['PCAs']
                new_descriptors = pd.np.empty((data.index.max()+1,len(particles),N_PCA))
                new_descriptors[:,:2,:] = PCAs[0].transform(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,N_PCA)
                new_descriptors[:,2:,:] = PCAs[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,N_PCA)
                descriptors = new_descriptors
                nb_features = N_PCA

    return metadata,data.join(pd.DataFrame({'descriptor':list(descriptors)}))

        
def create_data_NN(data,metadata):
    particles, scaler, test_size, H_HNN,H_ONN,O_HNN,O_ONN,N_NN = [metadata[x] for x in ['particles','scaler','test_size','H_HNN','H_ONN','O_HNN','O_ONN','N_NN']]
    local_dist_mat = np.empty((data.index.max()+1,len(particles),N_NN))
    nb_features = N_NN
    for i_time in tqdm.tqdm(range(data.index.max()+1)):
        dist_mat = cdist(data['molec'][i_time].get_positions(),data['molec'][i_time].get_positions())
        for i_part in range(len(particles)):
            if i_part in [0,1]:
                order = np.argsort(dist_mat[i_part])[1:]
                first_Hs = order[np.isin(order,range(2,7))][:O_HNN]
                first_Os = order[np.isin(order,range(0,2))][:O_ONN]
                local_dist_mat[i_time,i_part] = dist_mat[i_part,np.concatenate([first_Hs,first_Os])]
            else:
                order = np.argsort(dist_mat[i_part])[1:]
                first_Hs = order[np.isin(order,range(2,7))][:H_HNN]
                first_Os = order[np.isin(order,range(0,2))][:H_ONN]
                local_dist_mat[i_time,i_part] = dist_mat[i_part,np.concatenate([first_Hs,first_Os])]
                
                
    try:
        data['is_train']
    except KeyError:
        data['is_train'] = create_is_train(data.index.max()+1)
    else:
        pass

                
    if scaler == False:
        pass
        
    elif type(scaler) == type(None):
        local_dist_mat,scaler = scale_descriptors(data,local_dist_mat)

    else:
        local_dist_mat[:,0:2,:] = scaler[0].transform(local_dist_mat[:,0:2,:].reshape(local_dist_mat[:,0:2,:].shape[0]*2,nb_features)).reshape(local_dist_mat.shape[0],2,nb_features)
        local_dist_mat[:,2:,:] = scaler[1].transform(local_dist_mat[:,2:,:].reshape(local_dist_mat[:,2:,:].shape[0]*5,nb_features)).reshape(local_dist_mat.shape[0],5,nb_features)
    
    metadata['scaler'] = scaler
    return data.join(pd.DataFrame({'descriptor':list(local_dist_mat)})),metadata


def select_best_params(descriptors,nb_features,N_PCA):
    
    pca_O = PCA(n_components=N_PCA).fit(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features))
    pca_H = PCA(n_components=N_PCA).fit(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features))
    precision = [np.cumsum(pca_O.explained_variance_ratio_)[-1],np.cumsum(pca_H.explained_variance_ratio_)[-1]]
    
    print("Precision of new features [O,H] = ",precision)
        
    return [pca_O, pca_H]
    
    
    
    
def create_data_ACSF(data, metadata):
    particles, scaler, test_size,rcut, nmax, lmax,N_PCA,sigma_SOAP = [metadata[x] for x in ['particles','scaler','test_size','rcut','nmax','lmax','N_PCA','sigma_SOAP']]
                     
                     
    acsf = ACSF(
    species=["H", "O"],
    rcut=9.0,
    g2_params=[[1, 0], [0.1, 0], [0.01, 0],[0.01, 0],[0.001, 0]],
    g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1],[0.1, 1, 1], [0.1, 2, 1], [0.1, 1, -1], [0.1, 2, -1],    [0.01, 1, 1], [0.01, 2, 1], [0.01, 1, -1], [0.01, 2, -1]]
    )


    nb_features = acsf.get_number_of_features()
    
    descriptors = pd.np.empty((data.index.max()+1,len(particles),nb_features))
    
    for i_time in tqdm.tqdm(range(data.index.max()+1)):
        descriptors[i_time] = acsf.create(data['molec'][i_time],positions=np.arange(len(particles)))
    
    #create training set
    try:
        data['is_train']
    except KeyError:
        data['is_train'] = create_is_train(data.index.max()+1)
    else:
        pass
    #selecting best params
    if N_PCA:
            try:
                metadata['PCAs']
            except KeyError:
                PCAs = select_best_params(descriptors[data['is_train'].values],nb_features,N_PCA) 
                new_descriptors = pd.np.empty((data.index.max()+1,len(particles),N_PCA))
                new_descriptors[:,:2,:] = PCAs[0].transform(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,N_PCA)
                new_descriptors[:,2:,:] = PCAs[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,N_PCA)
                descriptors = new_descriptors
                metadata['old_N_feature'] = nb_features
                nb_features = N_PCA
                metadata['PCAs'] = PCAs
                

            else:
                PCAs = metadata['PCAs']
                new_descriptors = pd.np.empty((data.index.max()+1,len(particles),N_PCA))
                new_descriptors[:,:2,:] = PCAs[0].transform(descriptors[:,:2,:].reshape(descriptors[:,:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,N_PCA)
                new_descriptors[:,2:,:] = PCAs[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,N_PCA)
                descriptors = new_descriptors
                nb_features = N_PCA

    else:
        pass 
    #scaling
    if scaler == False:
        pass
        
    elif type(scaler) == type(None):
        descriptors,scaler = scale_descriptors(data,descriptors)

    else:
        descriptors[:,0:2,:] = scaler[0].transform(descriptors[:,0:2,:].reshape(descriptors[:,0:2,:].shape[0]*2,nb_features)).reshape(descriptors.shape[0],2,nb_features)
        descriptors[:,2:,:] = scaler[1].transform(descriptors[:,2:,:].reshape(descriptors[:,2:,:].shape[0]*5,nb_features)).reshape(descriptors.shape[0],5,nb_features)

    metadata['scaler'] = scaler
    return data.join(pd.DataFrame({'descriptor':list(descriptors)})),metadata
    
    
    
    
    
    