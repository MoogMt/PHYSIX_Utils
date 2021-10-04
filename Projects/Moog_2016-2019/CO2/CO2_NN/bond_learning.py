# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 12:19:24 2021

@author: moogm
"""

# Load necessary libraries
#--------------------------------------
# Pandas -> Data handling
import pandas as pd 
# Numpy -> Vector/array operations
import numpy as np
# Sklearn -> Machine Learning/Preprocessing
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn import svm
from sklearn.model_selection import train_test_split
# Matplotlib -> Data Visualization
import matplotlib.pyplot as plt
#--------------------------------------

# Data pre-treatment
#============================================================================
# Read CSV files
#--------------------------------------------------------------------
# File path to the datasets
density_filepath="F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\8.82\\Trajectory_2\\input_dens.dat"
elf_filepath="F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\8.82\\Trajectory_2\\input_elf.dat"
# Reading data files 
df_density = pd.read_csv( density_filepath, sep=" " )
df_elf = pd.read_csv( elf_filepath, sep=" " )
#--------------------------------------------------------------------

# Getting data dimensions (same for both datasets)
#--------------------------------------------------------------------
nb_data_point = df_density.shape[0]
nb_cols       = df_density.shape[1]
#--------------------------------------------------------------------

# Set the names of the columns
#--------------------------------------------------------------------
# - First columns indicates the frame
names_meta=["Step","Carbon","Neighbor","Oxygen","C-O Distance"]
# - Create name for ELF/Density values on line
nb_meta_cols=len(names_meta)
names = names_meta
for i in range(nb_meta_cols+1,nb_cols):
    names.append("Line " + str(i-nb_meta_cols) )
# Add output name
names.append("Status")
# - Replace names of the columns in the DF
df_density = df_density.set_axis( names, axis=1, inplace=False )
df_elf     = df_elf.set_axis(     names, axis=1, inplace=False )
#--------------------------------------------------------------------

# Separating between data between trainable and not trainable
#--------------------------------------------------------------------
mask_data = df_density["Status"] != 0
# Separate between trainable and application sets
df_density_trainable = df_density[  mask_data ]
df_density_test      = df_density[ ~mask_data ]
df_elf_trainable = df_elf[  mask_data ]
df_elf_test      = df_elf[ ~mask_data ]
#--------------------------------------------------------------------

# Separate between bonded and not bonded sets
#--------------------------------------------------------------------
df_bonded_density     = df_density[ df_density["Prediction"] == 1 ]
df_non_bonded_density = df_density[ df_density["Prediction"] == -1 ]
df_bonded_elf     = df_elf[ df_elf["Prediction"] ==  1 ]
df_non_bonded_elf = df_elf[ df_elf["Prediction"] == -1 ]
#--------------------------------------------------------------------

# Separating input and output 
#--------------------------------------------------------------------
# - input 
density_X = df_density_trainable.iloc[:,nb_meta_cols:nb_cols-1]
elf_X     = df_elf_trainable.iloc[:,nb_meta_cols:nb_cols-1]
# - output 
density_y = df_density_trainable.iloc[:,nb_cols-1]
elf_y     = df_elf_trainable.iloc[:,nb_cols-1]
#==============================================================================

# Preprocessing
#==============================================================================
# Arguments
#--------------------
test_size = 0.4
random_state = 42
shuffle = True
#--------------------

# Functions
#--------------------
def fitMinMaxScaler( X ):
    # Argument 
    # - X: input data
    # Output 
    # - Fitted Scaler (MinMax)
    return preprocessing.MinMaxScaler().fit(X)
def fitStdScaler( X ):
    # Argument
    # - X: input data
    # Output
    # - Fitter Scaler
    return preprocessing.StandardScaler().fit(X)
#--------------------

# Content
#--------------------
# Separating Train and Test split
# - density
density_X_train, density_X_test, density_y_train, density_y_test = train_test_split( density_X, density_y, test_size=test_size, random_state=random_state, shuffle=shuffle)
# - elf
elf_X_train, elf_X_test, elf_y_train, elf_y_test = train_test_split( elf_X, elf_y, test_size=test_size, random_state=random_state, shuffle=shuffle)
#============================================================================


# Machine Learning Application
#============================================================================
# Functions
#----------------------------------------------------
    # Arguments:
    # - X: input data
    # - n_components: number of components to keep
    
    # Returns PCA object, fitted on data
    return PCA( n_components=n_components ).fit(X
def learnSVC( X, y, reg_param=1.0, kernel="rbf" ):
    # Arguments:
    # - X: input data
    # - y: output data
    # - kernel: type of the kernel amongst ‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’
    # Output
    # - SVC predictor trained
    
    return svm.SVC( C=reg_param, kernel=kernel).fit( X, y )
#----------------------------------------------------

# Content
#----------------------------------------------------
# Create ML operator
clf_dens = learnSVC( density_X_train, density_y_train, reg_param=1.5, kernel="linear" )
clf_elf  = learnSVC( elf_X_train, elf_y_test, reg_param=1.5, kernel="linear" )

# Predicting output for the whole set (train+test+application)
df_density["Prediction"] = clf_dens.predict( density_X )
df_elf["Prediction"] = clf_elf.predict( elf_X )

# Predicting on the train and test sets
# - Density
y_pred_train_dens = clf_dens.predict( density_X_train )
y_pred_test_dens  = clf_dens.predict( density_X_test  )
# - ELF
y_pred_train_elf = clf_elf.predict( elf_X_train )
y_pred_test_elf  = clf_elf.predict( elf_X_test  )
#============================================================================

# Post-Processing
#============================================================================
# Mask containing only bonded distances, for the extrapolating set
mask_bonded_dens = y_test_dens == 1
mask_bonded_elf  = y_test_elf  == 1
#============================================================================

# Visualization
#============================================================================

# Parameters
#--------------------------
nb_bins=50          # Number of bins for the histogram in the whole dataset
nb_bins_extrapol=20 # Number of bins for the histogram in the extrapolation dataset
#--------------------------

# Plots of the whole distance distributions
#-------------------------------------------
plt.figure(1)
plt.hist( [df_bonded_density["C-O Distance"], df_non_bonded_density["C-O Distance"] ], bins=nb_bins, density=True, histtype='barstacked' )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of C-O distances")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(2)
plt.hist( [df_bonded_elf["C-O Distance"], df_non_bonded_elf["C-O Distance"]], bins=nb_bins, density=True, histtype='barstacked' )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of C-O distances")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(8)
plt.hist( [df_bonded_density["C-O Distance"], df_non_bonded_density["C-O Distance"] ], bins=nb_bins, density=True )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of C-O distances")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(9)
plt.hist( [df_bonded_elf["C-O Distance"], df_non_bonded_elf["C-O Distance"]], bins=nb_bins, density=False )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of C-O distances")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(3)
plt.hist( [ df_bonded_density["C-O Distance"], df_bonded_elf["C-O Distance"]] , bins=nb_bins, density=True )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of distances\n between bonded C-O")
plt.legend(["Density","ELF"])
plt.show()
#-------------------------------------------
plt.figure(4)
plt.hist( [ df_non_bonded_density["C-O Distance"], df_non_bonded_elf["C-O Distance"]] , bins=nb_bins, density=True )
plt.xlabel("Distance (A)")
plt.ylabel("Distribution of distances\n between non-bonded C-O")
plt.legend(["Density","ELF"])
plt.show()
#-------------------------------------------


# Plots only on the extrapolating set
#-------------------------------------------
plt.figure(4)
plt.hist( [ df_elf_test[ mask_bonded_elf ]["C-O Distance"], df_elf_test[ ~mask_bonded_elf ]["C-O Distance"] ], bins=nb_bins_extrapol, density=True, histtype='barstacked' )
plt.ylabel("Distribution")
plt.xlabel("Distance (A)")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(5)
plt.hist( [ df_density_test[ mask_bonded_dens ]["C-O Distance"], df_density_test[ ~mask_bonded_dens ]["C-O Distance"] ], bins=nb_bins_extrapol, density=True, histtype='barstacked' )
plt.ylabel("Distribution")
plt.xlabel("Distance (A)")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(6)
plt.hist( [ df_elf_test[ mask_bonded_elf ]["C-O Distance"], df_elf_test[ ~mask_bonded_elf ]["C-O Distance"] ], bins=nb_bins_extrapol, density=True )
plt.ylabel("Distribution")
plt.xlabel("Distance (A)")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#-------------------------------------------
plt.figure(7)
plt.hist( [ df_density_test[ mask_bonded_dens ]["C-O Distance"], df_density_test[ ~mask_bonded_dens ]["C-O Distance"] ], bins=nb_bins_extrapol, density=True )
plt.ylabel("Distribution")
plt.xlabel("Distance (A)")
plt.legend(["Bonded","Non-Bonded"])
plt.show()
#============================================================================