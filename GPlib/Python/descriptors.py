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
    

    
    
    
    
    