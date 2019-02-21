GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=5001
d_lag=5
unit=0.005

V=8.82
T=3000

# for V in Volumes
#     for T in Temperatures

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")
#folder_out=string(folder_in)


print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
states, percent, state_matrix, type_states = assignDataToStates( data, size(types)[1], type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"final_state_matrix.dat"), state_matrix )


cut_off_states = 0.1
states, type_states = isolateSignificantStates( states, percent, cut_off_states, type_states )
state_matrix, percent = assignDataToStates( data , states , size(types)[1], type_states, type_atoms, false)
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

# end
# end
