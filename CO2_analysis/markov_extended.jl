GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion

# Folder for data


folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=1001
d_lag=2
unit=0.005

# for V in Volume
#     for T in Temperatures


T=3000
V=8.82

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

states, state_matrix, count_ = assignDataToStates( data, nb_types, type_atoms )
#writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
#writeStateMatrix( string(folder_out,"initial_state_matrix.dat"), state_matrix )



percent_ = count_[1]/sum(count_[1])*100

check=zeros(size(states[1])[1])
for i=1:size(check)[1]
    if percent_[i] > 1
        check[i] = 1
    end
end

file_out=open(string(folder_out,"states_saved.dat"),"w")
sum_=0
for i=1:size(states[1])[1]
    if check[i] > 0
        for j=1:size(states[1])[2]
            write(file_out,string(states[1][i,j]," "))
        end
        write(file_out,string(percent_[i],"\n"))
        global sum_ += percent_[i]
    end
end
close(file_out)

transitions_matrix=[transitionMatrix(states[1],state_matrix[1],nb_types,type_atoms,min_lag,max_lag,d_lag)]
transitions_matrix_CK=[chappmanKormologov( transitions_matrix[1])]

file_out=open(string(folder_out,"transition_saved.dat"),"w")
file_out_CK=open(string(folder_out,"transition_saved_ck.dat"),"w")
nb_states=size(transitions_matrix[1])[1]
nb_time=size(transitions_matrix[1])[3]
for tau=1:nb_time
    write(file_out,string(tau*unit*d_lag," "))
    write(file_out_CK,string(tau*2*unit*d_lag," "))
for i=1:size(transitions_matrix[1])[1]
    #if check[i] > 0
        for j=1:size(transitions_matrix[1])[2]
            #if check[j] > 0
                write(file_out,string(transitions_matrix[1][i,j,tau]," "))
                if tau < nb_time/2
                    write(file_out_CK,string(transitions_matrix_CK[1][i,j,tau]," "))
                end
            #end
        end

    #end
end
write(file_out,string("\n"))
write(file_out_CK,string("\n"))
end
close(file_out)
close(file_out_CK)


stateO=[3,3,0,0,0,0,0,0,0,0]
stateC=[0,0,0,0,0,2,2,0,0,0]

O_nb=0
for i=1:size(states)[1]
    dO=0
    for j=1:size(states)[2]
        dO += (states[i,j]-stateO[j])^2
    end
    if dO == 0
        O_nb = i
    end
end

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

#     end
# end
