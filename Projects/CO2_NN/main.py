#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:50:29 2019

@author: moogmt
"""

from import_data import import_data_and_after
import descriptors
from predictors import energy_predictor, energy_predictor_dropout
from graphics import make_report

import sys
import os
import numpy as np
import pandas as pd

import metadata as mtd


metadata=mtd.buildMetaData(2,3,4,5,6)

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

