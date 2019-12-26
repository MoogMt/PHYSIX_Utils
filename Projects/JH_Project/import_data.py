#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:47:20 2019

@author: julienh

Contains functions for importing datas
"""

import numpy as np
import tqdm 
import pandas as pd
from ase import Atoms
import scipy as sp


def import_data_and_after(metadata):
    directory,tot_time,particles,size_file,replace,LJ_pot = [metadata[x] for x in ['directory_to_input_data','total_time','particles','time_of_file','replace','LJ_pot']]
    
    N_part = len(particles)
    all_positions = np.empty((size_file,N_part,3))
    all_energies = np.empty(size_file)

    with open(directory+'zundel.xyz', "r") as file_pos, open(directory+'energy.out',"r") as file_energies :

        file_energies.readline()

        for i_time in tqdm.tqdm(range(size_file)):
            file_pos.readline()
            file_pos.readline()
            for i_atom in range(N_part):
                all_positions[i_time,i_atom] = file_pos.readline().split()[1:]
                
            all_energies[i_time] = file_energies.readline().split()[4]
    
    molecs = np.empty(tot_time,dtype=Atoms)
    chosen_times = np.random.choice(size_file,size=tot_time,replace=replace)
    positions = all_positions[chosen_times]
    energies = all_energies[chosen_times]
    for i_time in tqdm.tqdm(range(tot_time)):
        molecs[i_time] = Atoms(numbers=particles, positions=positions[i_time])   
        
    if LJ_pot:
        energies = dummy_energy(positions)/627.509391

    
    data = pd.DataFrame({"energy":energies*627.509391,'molec':molecs})
    return data

def import_data_behler(metadata):
    directory,tot_time,particles,size_file,replace = [metadata[x] for x in ['directory_to_input_data','total_time','particles','time_of_file','replace']]
    
    
def dummy_energy(positions):
    energy_LJ = np.empty(positions.shape[0])
    for i_time in tqdm.tqdm(range(energy_LJ.size)):
        all_dist = sp.spatial.distance.pdist(positions[i_time])
        energy_LJ[i_time] =LJ_energy(all_dist)
    return energy_LJ
        
def LJ_energy(all_dist):
    d=2.2
    return ((d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist) - (d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)*(d/all_dist)).sum()/5
    
