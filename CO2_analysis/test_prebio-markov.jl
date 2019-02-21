GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Cut-off distance for bonds
cut_off_bond = 1.75

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
V=16.36074

# Defining base folder
folder_in=string("/media/moogmt/Stock/Theo_equilibration/")

# File of trajectory to be read
file=string(folder_in,"01_equilibration.xyz")

# Folder where to put the out data
folder_out=string(folder_in,"Data/")

# Reading xyz
print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

# Computing neighbours
data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

# Creating states and assigning data
states, percent, state_matrix, type_states = assignDataToStates( data, nb_types, type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"initial_state_matrix.dat"), state_matrix )

# Isolating only pertinent states for transition study
states, type_states = isolateSignificantStates( states, percent, cut_off_states, type_states )
state_matrix, percent = assignDataToStates( data , states , nb_types , type_states, type_atoms, false)
writeStates(string(folder_out,"markov_final_states-",cut_off_states,".dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"final_state_matrix.dat"), state_matrix )

# Transition matrix study
transitions_matrix=transitionMatrix(states,state_matrix,type_states,nb_types,type_atoms,min_lag,max_lag,d_lag)
# test of Chappman Kolmogorov test (Markov Testing not necessary if one just wants to study the transitions)
transitions_matrix_CK=[chappmanKormologov( transitions_matrix[1])]
for i=2:nb_types
    push!( transitions_matrix_CK, chappmanKormologov( transitions_matrix[i] ) )
end

# Writting results
for i=1:nb_types
    writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[i],".dat"),transitions_matrix[i])
    writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[i],"-CK.dat"),transitions_matrix_CK[i])
end
