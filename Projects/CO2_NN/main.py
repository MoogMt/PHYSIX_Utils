#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:50:29 2019

@author: moogmt
"""

import numpy as np
import pandas as pd
import metadata as mtd
import datahand as dth

from ase.build import molecule
from dscribe.descriptors import SOAP
from ase.io import read

data_base  = "/media/moogmt/Elements/CO2/"

volume=8.82
temperature=3000 
run_nb=1

folder_in = data_base + str(volume) + "/" + str(temperature) + "K/" + str(run_nb) + "-run/"
folder_out = data_base + str(volume) + "/" + str(temperature) + "K/Data/"

# Reading trajectory
file_traj = folder_in + "TRAJEC.xyz"
traj = read(file_traj,index=':')
for i in range(len(traj)):
    traj[i].set_cell([volume, volume, volume])

n_atoms=np.shape(traj[0])[0]
metadata=mtd.buildMetaData("TRAJEC.xyz","ENERGIES",folder_out, temperature)

dth.readData(metadata)

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



