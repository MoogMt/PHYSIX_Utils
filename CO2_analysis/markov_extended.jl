GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
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

# for V in Volumes
#     for T in Temperatures

V=9.5
T=3000

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

# if ! isfile(file)
#     continue
# end

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

states, percent, state_matrix, type_states = assignDataToStates( data, nb_types, type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"initial_state_matrix.dat"), state_matrix )

stateO=[3,3,0,0,0,0,0,0]
stateC=[0,0,0,0,2,2,0,0]

O_nb=0
C_nb=0
for i=1:size(states)[1]
    dC=0
    dO=0
    for j=1:size(states[1])[1]
        dC += (states[i,j]-stateC)^2
        dO += (states[i,j]-stateO)^2
    end
    if dC == 0
        C_nb = i
    end
    if dO == 0
        O_nb = i
    end
end

# cut_off_states = 0.1
# states, type_states = isolateSignificantStates( states, percent, cut_off_states, type_states )
# state_matrix, percent = assignDataToStates( data , states , nb_types, type_states, type_atoms, false)
# writeStates(string(folder_out,"markov_final_states-",cut_off_states,".dat"),states,percent,types,type_states)
# writeStateMatrix( string(folder_out,"final_state_matrix.dat"), state_matrix )
#
# # Transition matrix study
# transitions_matrix=transitionMatrix(states,state_matrix,type_states,nb_types,type_atoms,min_lag,max_lag,d_lag)
# # test of Chappman Kolmogorov test (Markov Testing not necessary if one just wants to study the transitions)
# transitions_matrix_CK=[chappmanKormologov( transitions_matrix[1])]
# for i=2:nb_types
#         push!( transitions_matrix_CK, chappmanKormologov( transitions_matrix[i] ) )
# end
#
# # Writting results
# for i=1:nb_types
#         writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[i],".dat"),transitions_matrix[i])
#         writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[i],"-CK.dat"),transitions_matrix_CK[i])
# end

taumax=100
dtau=1
file_out=open(string(folder_out,"Test.dat"),"w")
for step=1:size(state_matrix)[2]
    print("Progress: ",step/size(state_matrix)[2]*100,"%\n")
    count=0
    for oxygen=1:nbC
        if state_matrix[nbC+oxygen,step] == O_nb
            for oxygen2=oxygen+1:nbO
                if state_matrix[nbC+oxygen2,step] == O_nb
                    for carbon=1:nbC
                        if cell_mod.distance(traj[step],cell,carbon,nbC+oxygen) < cut_off_bond && cell_mod.distance(traj[step],cell,carbon,nbC+oxygen2) < cut_off_bond
                            for carbon2=carbon+1:nbC
                                if cell_mod.distance(traj[step],cell,carbon2,nbC+oxygen) < cut_off_bond && cell_mod.distance(traj[step],cell,carbon2,nbC+oxygen2) < cut_off_bond
                                    count += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    write(file_out,string(step," ",count,"\n"))
end
close(file_out)
#
#     end
# end
