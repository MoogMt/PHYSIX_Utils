# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 12:19:24 2021

@author: moogm
"""

import pandas as pd
import numpy as np

density_filepath="F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\8.82\\Trajectory_2\\input_dens.dat"
elf_filepath="F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\8.82\\Trajectory_2\\input_elf.dat"


col=27
name_col=3

names = ["Structure step", "Carbon","Neighbor Nb" ]
for i in range(col-name_col):
    names.append(str(i))

df_density = pd.read_csv( density_filepath, sep=" ", names=names )
df_elf = pd.read_csv( elf_filepath, sep=" ", names=names )