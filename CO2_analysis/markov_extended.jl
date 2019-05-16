GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[9.35,9.375,9.4,9.5,9.8,10.0]
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


T=3000
V=9.8

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

states, state_matrix, count = assignDataToStates( data, nb_types, type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)
writeStateMatrix( string(folder_out,"initial_state_matrix.dat"), state_matrix )

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
