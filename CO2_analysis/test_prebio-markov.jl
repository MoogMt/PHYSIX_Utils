GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
#folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=2001
d_lag=2
unit=0.005

folder_in=folder_base
file=string(folder_base,"01_equilibration.xyz")

folder_out=folder_base
#folder_out=string(folder_in)

V=16.36074

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
states, percent, state_matrix, type_states = assignDataToStates( data, size(types)[1], type_atoms )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent,types,type_states)

cut_off_states = 10
states, type_states = isolateSignificantStates( states, percent,  type_states, cut_off_states )
state_matrix, percent, unused_percent = assignDataToStates( data , states , false)
writeStates(string(folder_out,"markov_final_states-",cut_off_states,".dat"),states,percent,types,type_list)

# Checking chappmanKormologov
transition_matrix = transitionMatrix( states, state_matrix, min_lag, max_lag, d_lag )
transition_matrix_CK = chappmanKormologov( transition_matrix )

nb_states=size(states)[1]

for j=1:nb_states
    file_out=open(string(folder_out,"markov_CK_test-All-",cut_off_bond,"-",j,"-part1.dat"),"w")
    for i=1:2:size(transition_matrix)[3]
        write(file_out,string(i*unit*d_lag," "))
        for k=1:nb_states
            write(file_out,string(transition_matrix[j,k,i]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
end
for j=1:nb_states
    file_out=open(string(folder_out,"markov_CK_test-All-",cut_off_bond,"-",j,"-part2.dat"),"w")
    for i=1:size(transition_matrix_CK)[3]
        write(file_out,string(2*i*unit*d_lag," "))
        for k=1:nb_states
            write(file_out,string(transition_matrix_CK[j,k,i]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
end
