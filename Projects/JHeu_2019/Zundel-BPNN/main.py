#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 15:13:39 2019

@author: julienh
"""

from import_data import import_data_and_after
import descriptors
from predictors import energy_predictor, energy_predictor_dropout
from graphics import make_report


import datetime
import sys
import os
import numpy as np
import pandas as pd

metadata ={
#General
'temperature':100,
'time_of_file':999_990,
'total_time':100_000,
'replace':False,
'particles':[8,8,1,1,1,1,1],   
'system_name':'Zundel',
'energy_calculation_method':'Coupled clusters',
'datetime':datetime.datetime.now().isoformat(),

 #Descriptor
'descriptor_type':'NN',               #Choice of descriptor
'N_PCA':False,                            #False if no PCA, a number if PCA to select N first axis
#SOAP
'scaler':None,             #False if no scaling, None if scaling
'sigma_SOAP':0.8,
'rcut':9.0, 
'nmax':3,
'lmax':2,
#Nearest neighbors/PIV:
'N_NN':4,
'H_HNN':2,
'H_ONN':2,
'O_HNN':3,
'O_ONN':1,

#Neural Net
'test_size':0.10,
'activation_function':'tanh',
'loss_function':'mean_squared_error',
'optimizer':'adam', 
'epochs':1000,
'patience':100,
'N_nodes_H':[80,80], 
'N_nodes_O':[80,80], 


#Paths and imports
'test_data':False,
'LJ_pot':False,

"import_from":None,         #To import from precedent NN, give datetime, or None
"path_to_output":"/home/julienh/Desktop/outputNN/",
'folder_to_downfolded':'/home/julienh/Desktop/Julien/',  

}
metadata['directory_to_input_data'] = '/home/julienh/Desktop/data/coupled_cluster/'+str(metadata['temperature'])+'K/'
#metadata['directory_to_input_data'] = '/home/julienh/Desktop/data/training-data_H2O/input.data.BLYP'

os.mkdir(metadata['path_to_output']+metadata['datetime'])


#Importing data
if metadata['test_data']:
    data = pd.read_pickle("/home/julienh/Desktop/data/coupled_cluster/10_000_at_100K.pkl")
else:
    data = import_data_and_after(metadata)
#Choice of descriptor
data, metadata = getattr(descriptors, 'create_data_'+metadata['descriptor_type'])(data,metadata)

metadata["N_feature"]=data["descriptor"][0].shape[1]

#data, metadata["model"], metadata['mean_error'] = energy_predictor(data,metadata)
data, metadata["model"], metadata['mean_error'] = energy_predictor_dropout(data,metadata)


make_report(data, metadata)

np.save(metadata['path_to_output']+metadata['datetime']+'/metadata',metadata)


#Distinguer H_middle pour le scaling => COMMENT FAIRE LE SCALING ???

#Prevenir Ruben que son truc prend pas en compte qui est le h middle

#Parler à matthieu de comment scale ses data...

#Kalmann filters ?

#TESTER AVEC ET SANS AX+
#TESTER AVEC POIDS A LA FIN ET NON SUM



#DIRE EXP











#angular : all  same cut off 6ang
#1 2 4 16 
#1 2 for lambda


# largest extension // shortest distance in the dataset //
#

