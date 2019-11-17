GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.6]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]


nbC=32
nbO=nbC*2

cut_off_bond = 1.6
max_neigh=5

min_lag=1
max_lag=5001
d_lag=5
unit=0.005

V=8.82
T=3000

# for T in Temperatures
#
#     file_coordinanceC=open(string(folder_base,"coordinancesC-",T,".dat"),"w")
#     file_coordinanceO=open(string(folder_base,"coordinancesO-",T,".dat"),"w")
#     file_CO2=open(string(folder_base,"CO2-",T,".dat"),"w")
#     file_charged=open(string(folder_base,"charged-",T,".dat"),"w")
#     file_exchange=open(string(folder_base,"exchange-",T,".dat"),"w")
#     file_exchange_count=open(string(folder_base,"exchange_count-",T,".dat"),"w")
#
#     for V in Volumes

print("V=",V," T=",T,"\n")

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

# if ! isfile(file)
#     continue
# end


# fileP="Avg_Pressure-BootStrap-nboot_1000.dat"
# if ! isfile( string(folder_out,fileP) )
#     continue
# end

# file_p=open(string(folder_out,fileP))
# lines=readlines(file_p);
# close(file_p)
# P=parse(Float64,split(lines[1])[2])

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]
states, state_matrices, counts = assignDataToStates( data, nb_types, type_atoms )

min_lag=1
max_lag=5001
d_lag=2
unit=0.005

transition_matrix = transitionMatrix(states[1],state_matrices[1],nb_types,type_atoms,min_lag,max_lag,d_lag)
file_out=open(string(folder_out,"TM-",cut_off_bond,"-",d_lag,"-",max_lag,".dat"),"w")
for time_=1:size(transition_matrix)[3]
    write(file_out,string(time_*unit*d_lag," ",))
    for state_origin=1:size(transition_matrix)[1]
        for state_dest=1:size(transition_matrix)[2]
            write(file_out,string(transition_matrix[state_origin,state_dest,time_]*100," "))
        end
    end
    write(file_out,string("\n"))
end
close(file_out)




lifetime=0
lag=10
state=1
for time_=1:size(transition_matrix)[3]
    global lifetime += 0.5*unit*d_lag*time_*lag*(transition_matrix[state,state,lag]^time_)
end


# case_states=size(states[1])[2]
# nb_states=size(states[1])[1]
# normC=4
# normO=2
# oks_states=zeros(nb_states)
# for state=1:nb_states
#     part1=4
#     part2=0
#     for i=1:case_states
#         if states[1][state,i] > 0
#             part1 -= 1
#             if i > max_neigh # O
#                 part2 += (normO-states[1][state,i])
#             else # C
#                 part2 += (normC-states[1][state,i])
#             end
#         end
#     end
#     if part1 == part2
#         oks_states[state] = 1
#     end
# end
#
# nb_states_exc=0
# count_states_all=0
# for i=1:nb_states
#     if oks_states[i] != 1
#         nb_states_exc += percent[1][i]
#     end
#     count_states_all+= percent[1][i]
# end
# percent_excited=nb_states_exc/count_states_all
#
# write(file_charged,string(P," ",percent_excited,"\n"))

# coordinancesC=zeros(max_neigh)
# for state=1:nb_states
#     nb=0
#     for i=1:case_states
#         if states[1][state,i] > 0
#             nb += 1
#         end
#     end
#     coordinancesC[nb] += percent[1][state]
# end
# for i=1:max_neigh
#     coordinancesC[i] /= count_states_all
# end
# coordinancesC *= 100
#
# write(file_coordinanceC,string(P," "))
# for i=1:max_neigh
#     write(file_coordinanceC,string(coordinancesC[i]," "))
# end
# write(file_coordinanceC,"\n")
#
# percent_CO2=0
# stateC=[0,0,0,0,0,1,1,0,0,0]
# for state=1:nb_states
#     found=true
#     for i=1:case_states
#         if stateC[i] != states[1][state,i]
#             found=false
#         end
#     end
#     if found
#         percent_CO2=percent[1][state]/count_states_all*100
#         break
#     end
# end
# write(file_CO2,string(P," ",percent_CO2,"\n"))
#
# stateC=[0,0,0,0,0,2,2,1,0,0]
# nb_state_dimer=0
# for state=1:nb_states
#     found=true
#     for i=1:case_states
#         if stateC[i] != states[1][state,i]
#             found=false
#         end
#     end
#     if found
#         global nb_state_dimer=state
#         break
#     end
# end
#
# stateO=[3,3,0,0,0,0,0,0,0,0]
# nb_state_dimerO=0
# for state=1:nb_states
#     found=true
#     for i=1:case_states
#         if stateO[i] != states[2][state,i]
#             found=false
#         end
#     end
#     if found
#         global nb_state_dimerO=state
#         break
#     end
# end
#

#
# dt=250
# time_stock=[]
# count_exchange=0
# carbon_stock=[]
# for carbon=1:nbC
#     record=zeros(nb_steps)
#     for step_=1:nb_steps
#         found_one=false
#         if state_matrix[1][carbon,step_] == nb_state_dimer
#             # Check the exchange in two step:
#             max_neigh2=zeros(2)
#             count_bond_main=0
#             for atom=1:nb_atoms
#                 if atom == carbon
#                     continue
#                 end
#                 step_before=step_-1
#                 if step_before <= 0
#                     step_before=1
#                 end
#                 step_after = step_+dt
#                 if step_after > nb_steps
#                     step_after=nb_steps
#                 end
#                 if cell_mod.distance(traj[step_before],cell,carbon,atom) < cut_off_bond && cell_mod.distance(traj[step_after],cell,carbon,atom) < cut_off_bond
#                     if traj[step_].names[atom] == "O"
#                         if state_matrix[2][atom-nbC,step_] == nb_state_dimerO
#                             count_bond_main += 1
#                         end
#                     end
#                 end
#             end
#             if count_bond_main > 0 && record[step_] < 1
#                 for i=1:dt
#                     if step_+dt < nb_steps
#                         record[step_+i] += 1
#                     else
#                         break
#                     end
#                 end
#                 global count_exchange += 1
#                 push!(time_stock,step_)
#                 push!(carbon_stock,carbon)
#                 found_one=true
#             end
#         end
#         if found_one
#             if step_+dt < nb_steps
#                 step_ += dt
#             else
#                 break
#             end
#         end
#     end
# end
#
#
# count_exchange2=0
# carbon_exchange2=zeros(nbC)
# time_stock2=[]
# carbon_stock2=[]
# check=zeros(size(carbon_stock)[1])
# for i=1:size(carbon_stock)[1]-1
#     print("size:",size(check)[1],"\n")
#     found=false
#     for j=i+1:size(carbon_stock)[1]
#         if check[j] == true
#             continue
#         end
#         if time_stock[i] == time_stock[j]
#             check[i] = 1
#             check[j] = 1
#             push!(time_stock2,time_stock[i])
#             push!(carbon_stock2,carbon_stock[i])
#             push!(time_stock2,time_stock[j])
#             push!(carbon_stock2,carbon_stock[j])
#             found=true
#             global count_exchange2 += 1
#             if carbon_exchange2[carbon_stock[j]] < 1
#                 carbon_exchange2[carbon_stock[j]] += 1
#             end
#             if carbon_exchange2[carbon_stock[i]] < 1
#                 carbon_exchange2[carbon_stock[i]] += 1
#             end
#             break
#         end
#     end
# end
# nb_carbons=sum(carbon_exchange2)
# frac_carbon=nb_carbons/nbC*100

#         write(file_exchange,string(P," ",frac_carbon,"\n"))
#         write(file_exchange_count,string(P," ",count_exchange2,"\n"))
#
#     end
#
#     close(file_coordinanceC)
#     close(file_coordinanceO)
#     close(file_CO2)
#     close(file_charged)
#     close(file_exchange)
#     close(file_exchange_count)
#
# end
