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

from dscribe.descriptors import SOAP
from dscribe.descriptors import ACSF

from scipy.spatial.distance import cdist

def createDescriptorsSOAP(data, metadata):   
    # Prepping SOAP descriptor structure
    soap = SOAP( species=np.unique(metadata['n_atoms']), sigma=metadata['sigma_SOAP'], periodic=metadata['periodic'], rcut=metadata['cutoff_SOAP'], nmax=metadata['nmax'], lmax=metadata['lmax'],sparse=metadata['sparse_SOAP'])
    metadata['n_features'] = soap.get_number_of_features()    
    # Computing descriptors
    descriptors = pd.np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_features']))
    for index in tqdm.tqdm(range(metadata['total_set_size'])):
        descriptors[index] = soap.create(data['structures'][index],positions=np.arange(metadata['n_atoms']))
    return metadata, data.join(pd.DataFrame({'descriptor':list(descriptors)}))

def createDescriptorsACSF(data, metadata):         
    # Prepping ACSF descriptor structure
    acsf = ACSF(species=metadata['species'],rcut=metadata['cutoff_acsf'],g2_params=metadata['g2_params'],g4_params=metadata['g3_params'])
    metadata['n_features'] = acsf.get_number_of_features()    
    # Computing descriptors
    descriptors = pd.np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_features']))    
    for i_time in tqdm.tqdm(range(metadata['total_set_size'])):
        descriptors[i_time] = acsf.create(data['molec'][i_time],positions=np.arange(metadata['n_atoms']))
    return metadata, data.join(pd.DataFrame({'descriptor':list(descriptors)}))
    
def create_data_NN(data,metadata):
    local_dist_mat = np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_nearest_neigh']))
    metadata['n_features'] = metadata['n_nearest_neigh']
    for index in tqdm.tqdm(range(metadata['total_set_size'])):
        dist_mat = cdist(data['structures'][index].get_positions(),data['structures'][index].get_positions())
        for atom in range(metadata['n_atoms']):
            if atom in [0,1]:
                order = np.argsort(dist_mat[atom])[1:]
                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
            else:
                order = np.argsort(dist_mat[atom])[1:]
                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
    return metadata, data.join(pd.DataFrame({'descriptor':list(local_dist_mat)}))
    
    
    
    
    