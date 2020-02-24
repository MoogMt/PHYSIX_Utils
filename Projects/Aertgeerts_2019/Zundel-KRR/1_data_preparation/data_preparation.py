#Python Imports

#ASE imports
import ase
from ase.build import molecule
from ase import Atoms
import pandas as pd
import os
import pickle

#Choose the system

#Import positions
#All columns must have names
#Potential energy column must be called 'potential_energy'
positions_path = '../0_raw_data/Zundel_CCFF_50K/zundel.xyz'

#Import energies
energies_path = '../0_raw_data/Zundel_CCFF_50K/energy.out'

#set metadata
metadata = {
    'position_unit': 'Bohr',
    'energy_unit': 'Kcal/mol',
    'energy_calculation_method' : 'Coupled Cluster',
    'temperature_of_simulation' : '50K',
    'system_name' : 'Zundel Ion',
    'timestep_between_data' : '20fs'
}

#Set destiantaion
save_dir = 'prepared_data'
filename = 'CCFF_zundel_50K_cleaned'


#Limit data range of the import
max_number_of_steps = 1000000
#Choose how much data should be kept
#iprint = 1 means all data is kept
#iprint = 20 one in twenty data entries are kept
iprint = 20
index_of_positions = slice(0,max_number_of_steps,iprint)

#Read the xyz file
positions = ase.io.read(positions_path, index=index_of_positions)
energies = pd.read_csv(energies_path, nrows=max_number_of_steps ,delim_whitespace = True)

#Select every iprint row
energies = energies[index_of_positions]
energies = energies.reset_index()

#Create the data file
data_imported = pd.DataFrame()
data_imported['configuration'] = positions
data_imported['energy'] = energies['potential_energy']

data = data_imported.copy()

#Only keep interesting columns
data = data[['configuration','energy']]

#Convert energy to Kcal/mol
data['energy'] = data['energy'].multiply(627.509)

#add the medatada
data['metadata'] = None
data['metadata'].iat[0] = metadata

#save the pickle file
import os
import pickle
save_path = os.path.join(save_dir,filename)
pickle.dump(data, open(save_path, 'wb'))

#Test if the dataframe was correctly saved
test = pickle.load(open(save_path, 'rb'))
test.head()

#Test metadata overview
test['metadata'][0]
