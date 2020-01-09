#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:49:53 2019

@author: julienh with modifications from moogmt
"""

import tqdm
import pandas as pd
import numpy as np

from dscribe.descriptors import SOAP
from dscribe.descriptors import ACSF

# - SOAP
#=============================================================================#

default_sigma_SOAP  = 0.05  # Sigma for the gaussian density
default_cutoff_SOAP = 1.001 # Angstroms
default_nmax_SOAP   = 1     # nmax for the radial expansion
default_lmax_SOAP   = 0     # lmax for the angular expansion
default_sparse_SOAP = False # Sparse, avoid 0

def checkSOAPParams( metadata, verbose ):
    if metadata['sigma_SOAP'] < 0:
        if verbose:
            print("Invalid value of the sigma for SOAP: sigma=",metadata['sigma_SOAP'],"\n")
        return False
    if metadata['cutoff_SOAP'] < 1.0 :
        if verbose:
            print("Invalue value of the cut_off for SOAP: cut off = ",metadata['cutoff_SOAP'],"\n")
        return False
    if metadata['nmax_SOAP'] < 1.0 :
        if verbose:
            print("Invalue value of the cut_off for SOAP: nmax = ",metadata['cutoff_SOAP'],"\n")
        return False
    if metadata['lmax_SOAP'] < 0. or metadata['lmax_SOAP'] > metadata['nmax_SOAP'] : 
        if verbose: 
            print("Invalue value of the cut_off for SOAP: lmax = ",metadata['lmax_SOAP'],"\n")
        return False
    return True

def createDescriptorsSOAP(data, metadata,
                          sigma_SOAP=default_sigma_SOAP, 
                          cutoff_SOAP=default_cutoff_SOAP, 
                          nmax_SOAP=default_nmax_SOAP, 
                          lmax_SOAP=default_lmax_SOAP,
                          sparse_SOAP=default_sparse_SOAP,
                          ):
    # Updating metadata
    metadata['sigma_SOAP']  = sigma_SOAP
    metadata['cutoff_SOAP'] = cutoff_SOAP
    metadata['nmax_SOAP']   = nmax_SOAP
    metadata['lmax_SOAP']   = lmax_SOAP
    metadata['sparse_SOAP'] = sparse_SOAP
    if not checkSOAPParams( metadata, metadata['verbose']): 
        return False, False
    # Prepping SOAP descriptor structure
    soap = SOAP( species=metadata['species'], sigma=metadata['sigma_SOAP'], periodic=metadata['periodic'], rcut=metadata['cutoff_SOAP'], nmax=metadata['nmax_SOAP'], lmax=metadata['lmax_SOAP'],sparse=metadata['sparse_SOAP'] )
    metadata['n_features'] = soap.get_number_of_features()
    # Computing descriptors
    descriptors = pd.np.empty((np.shape(data)[0],metadata['n_atoms'],metadata['n_features']))
    for index_structure in tqdm.tqdm(range(np.shape(data)[0])): # This whole thing is a bit insane, we could directly build it propertly for the training and avoid nonsentical stuff...
        for index_atom in range(metadata['n_atoms']):
            descriptors[index_structure,index_atom,:] = soap.create(data['structures'][index_structure],positions=[index_atom],)
    return metadata, descriptors
#=============================================================================#

# ACSF
# TO BE REDONE AT A LATER POINT
#=============================================================================#
default_cutoff_ACSF = 0.                          # Cut-off for the ACSF function
default_n_acsf = 0                                # Number of ACSF functions per atoms
default_g2_params = np.zeros((default_n_acsf,2))  # Parameters for the g2 functions for the ACSF functions
default_g3_params = np.zeros((default_n_acsf,3))  # Parameters for the g3 functions for the ACSF functions
def createDescriptorsACSF(data, metadata,
                          cutoff_ACSF=default_cutoff_ACSF,
                          g2_params=default_g2_params,
                          g3_params=default_g3_params,
                          ):   
    # Updating metadata      
    
    # Prepping ACSF descriptor structure
    acsf = ACSF(species=metadata['species'],rcut=metadata['cutoff_acsf'],g2_params=metadata['g2_params'],g4_params=metadata['g3_params'])
    metadata['n_features'] = acsf.get_number_of_features()    
    # Computing descriptors
    descriptors = pd.np.empty((metadata['train_set_size'],metadata['n_atoms'],metadata['n_features']))    
    for index_structure in tqdm.tqdm(range(metadata['train_set_size'])):
        for index_atom in tqdm.tqdm(range(metadata['n_atoms'])):
            descriptors[index_structure,index_atom,:] = acsf.create(data['structures'][index_structure],positions=[index_atom])
    return metadata, descriptors
#=============================================================================#    

#from scipy.spatial.distance import cdist

#def create_data_NN(data,metadata):
#    local_dist_mat = np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_nearest_neigh']))
#    metadata['n_features'] = metadata['n_nearest_neigh']
#    for index in tqdm.tqdm(range(metadata['total_set_size'])):
#        dist_mat = cdist(data['structures'][index].get_positions(),data['structures'][index].get_positions())
#        for atom in range(metadata['n_atoms']):
#            if atom in [0,1]:
#                order = np.argsort(dist_mat[atom])[1:]
#                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
#                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
#                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
#            else:
#                order = np.argsort(dist_mat[atom])[1:]
#                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
#                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
#                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
#    return metadata, data.join(pd.DataFrame({'descriptor':list(local_dist_mat)}))
#    
    
    
    
    