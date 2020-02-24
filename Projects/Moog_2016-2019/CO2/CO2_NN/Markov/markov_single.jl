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



# for V in Volume
#     for T in Temperatures


T=2500
V=9.4

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
percent_ = count_[1]/sum(count_[1])*100

min_lag=1
max_lag=1001
d_lag=2
unit=0.005

for carbon=1:nbC
    transitions_matrix=transitionMatrix(states[1],state_matrix[1][carbon,:],nb_types,type_atoms,min_lag,max_lag,d_lag)
    transitions_matrix_CK=chappmanKormologov( transitions_matrix)
    file_out=open(string(folder_out,"transsition_C_",carbon,".dat"),"w")
    for i=1:size(transitions_matrix)[3]
        Base.write(file_out,string(i*unit*d_lag," "))
         for j=1:size(transitions_matrix)[1]
             for k=1:size(transitions_matrix)[2]
                 Base.write(file_out,string(transitions_matrix[j,k,i]," "))
             end
        end
        Base.write(file_out,string("\n"))
    end
    close(file_out)
    file_out=open(string(folder_out,"transsition_C_",carbon,"_CK.dat"),"w")
    for i=1:size(transitions_matrix_CK)[3]
         Base.write(file_out,string(i*2*unit*d_lag," "))
         for j=1:size(transitions_matrix)[1]
             for k=1:size(transitions_matrix)[2]
                 Base.write(file_out,string(transitions_matrix_CK[j,k,i]," "))
             end
         end
         Base.write(file_out,string("\n"))
    end
    close(file_out)
end
