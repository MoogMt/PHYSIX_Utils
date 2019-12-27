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

def choseTrainDataByIndex(metadata,structures,energies,chosen_index):    
    structures = structures[chosen_index]
    energies = energies[chosen_index]
    if type(structures[0]) == ase.atoms.Atoms :
        structuresAtoms=np.empty(metadata['size_train_set'],dtype=ase.atoms.Atoms)
        for i in range(metadata['size_train_set']) :
            structuresAtoms[i] = ase.atoms.Atoms(numbers=metadata['n_atoms'], positions=structures[i,:] )   
        return pd.DataFrame({"energy":energies,'structure':structuresAtoms})
    return pd.DataFrame({"energy":energies,'structure':structures})


# 
def choseTrainDataRandom(metadata,structures,energies):    
    total_size_set=len(energies)
    chosen_index = np.random.choice(total_size_set,size=metadata['size_train_set'],replace=metadata['replace'])
    return choseTrainDataByIndex(metadata,structures,energies,chosen_index)

