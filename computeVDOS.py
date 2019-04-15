#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 21:52:23 2019

@author: moogmt
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

folder=str("/home/moogmt/CO2/CO2_AIMD/")

V=8.82
T=3000

# Recupere les donnees
f = open(str(folder+V+"/"+T+"K/"+"TRAJEC_wrapped.xyz"), "r")
f.close()