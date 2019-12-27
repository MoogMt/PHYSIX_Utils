#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 22:03:09 2019

@author: moogmt

Functions that deals with XYZ file 
Mostly built-in function in the ASE library 
# Remember to add proper ASE citations...
"""

import ase
import numpy as np


def getNbStep( file_path ):
    return len(ase.io.read(file_path),index=':')

def read( file_path ):
    return ase.io.read(file_path,index=':')

def readAsArray( file_path ):
    traj=ase.io.read(file_path,index=':')
    traj_ = np.zeros( len(traj),3)
    for i in range( len(traj)):
        traj_[i,:] = traj[0].positions
    return 
