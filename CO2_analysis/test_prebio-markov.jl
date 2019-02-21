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
d_lag=5
unit=0.005

# Volume of the cell (only orthorombic is implemented yet)
V=16.36074

# for V in Volumes
#     for T in Temperatures

folder_in=string("/media/moogmt/Stock/Theo_equilibration/")
file=string(folder_in,"01_equilibration.xyz")

folder_out=string(folder_in,"Data/")

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

states, percent, state_matrix, type_states = assignDataToStates( data, nb_types, type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"initial_state_matrix.dat"), state_matrix )

states, type_states = isolateSignificantStates( states, percent, cut_off_states, type_states )
state_matrix, percent = assignDataToStates( data , states , nb_types , type_states, type_atoms, false)
writeStates(string(folder_out,"markov_final_states-",cut_off_states,".dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"final_state_matrix.dat"), state_matrix )

transitions_matrix=transitionMatrix(states,state_matrix,type_states,nb_types,type_atoms,min_lag,max_lag,d_lag)
transitions_matrix_CK=[transitions_matrix[1]]
for i=1:nb_types
    chappmanKormologov( transition_matrix[1] )
end

# end
# end
