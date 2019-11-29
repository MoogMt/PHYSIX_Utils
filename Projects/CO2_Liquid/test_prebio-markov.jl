GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Cut-off distance for bonds
cut_off_bond = 2.2

# Maximum of number of neighbor for a given atom
max_neigh=5

# Cut-off to select states
cut_off_states = 0.1

# Parameters for the autocorrelation for the transition matrix
min_lag=1
max_lag=5001
d_lag=2
unit=0.001

# Volume of the cell (only orthorombic is implemented yet)
V=13.36074

# Defining base folder
folder_in=string("/media/moogmt/Stock/Theo/Equilibration/")

# File of trajectory to be read
file=string(folder_in,"01_equilibration.xyz")

# Folder where to put the out data
folder_out=string(folder_in,"Data/")

# Reading xyz
print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

folder_out=string(folder_in,"Data/")

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

states, state_matrices, counts_states = assignDataToStates( data, nb_types, type_atoms )

for type=1:2
writeStates( string(folder_out,"markov_intial_states-",types[type],"-",cut_off_bond,".dat"),states[type],counts_states[type],types)
writeStateMatrix( string(folder_out,"initial_state_matrix-",types[type],"-",cut_off_bond,".dat"), state_matrices[type] )
end

# Transition matrix study
transitions_matrix=[transitionMatrix(states[1],state_matrices[1],nb_types,type_atoms,min_lag,max_lag,d_lag)]
transitions_matrix_CK=[chappmanKormologov( transitions_matrix[1])]
for type=2:nb_types
push!( transitions_matrix, transitionMatrix( states[type], state_matrices[type], nb_types, type_atoms, min_lag, max_lag, d_lag) )
push!( transitions_matrix_CK, chappmanKormologov( transitions_matrix[type] ) )
end

for type=1:nb_types
writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[type],".dat"),transitions_matrix[type])
writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[type],"-CK.dat"),transitions_matrix_CK[type])
end
