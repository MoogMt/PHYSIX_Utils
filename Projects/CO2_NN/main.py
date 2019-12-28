#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:50:29 2019

@author: moogmt
"""

import numpy as np
import nnmetadata as mtd
import nndatahand as dth
import filexyz as xyz
import cpmd 

from dscribe.descriptors import SOAP

data_base  = "/media/moogmt/Elements/CO2/"

volume=8.82
temperature=3000 
run_nb=1

folder_in = data_base + str(volume) + "/" + str(temperature) + "K/" + str(run_nb) + "-run/"
folder_out = data_base + str(volume) + "/" + str(temperature) + "K/Data/"

file_traj = folder_in + "TRAJEC.xyz"
file_energies = folder_in + "ENERGIES"

metadata=mtd.buildMetaData(file_traj,file_energies,folder_out, temperature)
if not mtd.checkMetaDataIO(metadata,True):
    exit

# Reading trajectory
traj = xyz.readPbcCubic( file_traj, volume )
# Reading ENERGIES file
energies=cpmd.readPotEnergyCPMD( file_energies )

metadata['train_fraction'] = 0.2

data_train = dth.choseTrainDataRandom(metadata,traj,energies)

data, metadata = desc.createDescriptors(metadata)
network,metadata = nn.build_network(metadata)
network,metadata = nn.train_network(metadata)
output, metadata = nn.test_network(metadata)
io.write_output(metadata,output)

# Setting SOAP
species = ["C", "O"]
rcut = 3.0
nmax = 3
lmax = 3

# Setting up the SOAP descriptor
soap = SOAP(
    species=species,
    periodic=True,
    rcut=rcut,
    nmax=nmax,
    lmax=lmax,
)

soap_co2= soap.create(traj[0], positions=[[2.0,10,15]],n_jobs=8)



