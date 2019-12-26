#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 12:15:45 2019

@author: moogmt inspired by julienh

Contains functions for handling data
"""

import numpy as np
import tqdm 
import pandas as pd
from ase import Atoms

def readEnergiesCPMD(file_path):
    with open

def readData(metadata):
    directory,tot_time,particles,size_file,replace,LJ_pot = [metadata[x] for x in ['directory_to_input_data','total_time','particles','time_of_file','replace','LJ_pot']]

    traj = ase.io.read(file_input,index=':')
    
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
            
    data = pd.DataFrame({"energy":energies,'molec':molecs})
    return data
