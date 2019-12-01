GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"clustering.jl"))
include(string(CO2folder,"markovCO2.jl"))


# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"

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
#
# if ! isfile(file)
#     continue
# end
folder_out=string(folder_in,"Data/")

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)
nb_steps=size(traj)[1]

n_neighbor=4

filed=open(string(folder_out,"massive_data.dat"),"w")
for step=1:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")
    for carbon=1:nbC
        if state_matrices[1][carbon,step] == 1
            distances=zeros(nbO)
            for oxygen=1:nbO
                distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
            end
            index_sort=sortperm(distances)
            distances_2nd=ones(n_neighbor)*V
            for index_=1:n_neighbor
                for carbon2=1:nbC
                    if carbon2 != carbon
                        distanceOC=cell_mod.distance(traj[step],cell,carbon2,nbC+index_sort[index_])
                        if distances_2nd[index_] > distanceOC
                            distances_2nd[index_] = distanceOC
                        end
                    end
                end
            end
            write(filed,string(carbon," "))
            for k=1:4
                write(filed,string(distances[index_sort[k]]," ",distances_2nd[k]," "))
            end
            write(filed,string("\n"))
        end
    end
end
close(filed)



data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

states, state_matrices, counts_states = assignDataToStates( data, nb_types, type_atoms )

for type=1:2
    writeStates( string(folder_out,"markov_intial_states-",types[type],"-",cut_off_bond,".dat"),states[type],counts_states[type],types)
    writeStateMatrix( string(folder_out,"initial_state_matrix-",types[type],"-",cut_off_bond,".dat"), state_matrices[type] )
end


# Basic analysis of states


# Transition matrix study
transitions_matrix=[transitionMatrix(states[1],state_matrices[1],nb_types,type_atoms,min_lag,max_lag,d_lag)]
transitions_matrix_CK=[chappmanKormologov( transitions_matrix[1])]
for type=2:nb_types
    push!( transitions_matrix, transitionMatrix( states[type], state_matrices[type], nb_types, type_atoms, min_lag, max_lag, d_lag) )
    push!( transitions_matrix_CK, chappmanKormologov( transitions_matrix[type] ) )
end

using LinearAlgebra

function computeMFT( transition_matrix::Array{T1,3}, starting::T2, d_lag::T3, unit::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real }
    nb_states=size(transition_matrix)[1]
    nb_times=size(transition_matrix)[3]
    lifes=zeros( nb_states)

    Imatrix=Matrix{Float64}(I, 2, 2)
    W=transition_matrix[:,:,starting]^1000
    P=transition_matrix[:,:,starting]
    Z=inv(I-P+W)
    M=zeros(nb_states,nb_states)
    for i=1:nb_states
        for j=1:nb_states
            if i != j
                M[i,j]=(Z[j,j]-Z[i,j])/W[j,j]
            end
        end
    end

    return M
end

mean_first_passage_times1=computeMFT(transitions_matrix[1],1,d_lag,unit)


for type=1:nb_types
    writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[type],".dat"),transitions_matrix[type])
    writeTransitionsMatrix(string(folder_out,"TransitionsMatrix-",types[type],"-CK.dat"),transitions_matrix_CK[type])
end

# end
# end
# 2000 33.5ps
# 2500 9.51 ps
# 3000 3.9ps

# 8.82 3000 0.5ps
