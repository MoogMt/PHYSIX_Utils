#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 12:15:45 2019

@author:  julienh with modification from moogmt

Contains functions for handling data
"""

import numpy as np
import tqdm 
import pandas as pd
import ase

# Column code for ENERGIES file
cpmd_temperature_col=2
cpmd_pot_energy_col=3
cpmd_tot_energy_col=4
cpmd_msd_col=6
cpmd_scf_comptime=7

def getNbLineEnergiesCPMD(file_path):
    nb_line=0
    with open(file_path,"r") as f:
        f.readline()
        nb_line += 1
    return nb_line

def readPotEnergyCPMD(file_path):
    nb_point=getNbLineEnergiesCPMD(file_path)
    energies=np.zeros(nb_point)
    f=open(file_path,"r")
    for i in range (nb_point):
        energies[i] = f.readline().split()[cpmd_pot_energy_col] # Read Kohn-Sham energies (in Ry)
    f.close()
    return energies

def readEnergiesFile(file_path):
    nb_point=getNbLineEnergiesCPMD(file_path)
    data=np.zeros(nb_point,7)
    f=open(file_path,"r")
    for i in range (nb_point):
        data[i,:] = f.readline().split()[1:] # Read all data except first column (time)
    f.close()
    return data

def extractTemperature(data):
    return data[:,cpmd_tot_energy_col-1]

def extractPotentialEnergy(data):
    return data[:,cpmd_pot_energy_col-1]

def extractTotalEnergy(data):
    return data[:,cpmd_tot_energy_col-1]

def extractMSD(data):
    return data[:,cpmd_msd_col-1]

def extractSCFcomputationTime(data):
    return data[:,cpmd_scf_comptime-1]

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
