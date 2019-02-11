GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75

min_lag=1
max_lag=5001
d_lag=5
unit=0.005


T=3000
V=8.82


folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")
#folder_out=string(folder_in)

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_list=buildCoordinationMatrix( traj , cell , cut_off_bond )
states, percent, state_matrix = assignDataToStates( data, size(type)[1], type_list )
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent)

cut_off_states = 0.1
states = isolateSignificantStates( states, percent, cut_off_states )
state_matrix, percent, unused_percent = assignDataToStates( data , states , false)
transition_matrix = transitionMatrix( states, state_matrix, min_lag, max_lag, d_lag )
transition_matrix_CK = chappmanKormologov( transition_matrix )
writeStates(string(folder_out,"markov_final_states-",percent,".dat"),states,percent)

nb_states=size(states)[1]

for j=1:nb_states
    file_out=open(string(folder_out,"C_markov_CK_test-",cut_off_bond,"-",j,"-part1.dat"),"w")
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
    file_out=open(string(folder_out,"C_markov_CK_test-",cut_off_bond,"-",j,"-part2.dat"),"w")
    for i=1:size(transition_matrix_CK)[3]
        write(file_out,string(2*i*unit*d_lag," "))
        for k=1:nb_states
            write(file_out,string(transition_matrix_CK[j,k,i]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
end


# nb_states=size(transition_matrix)[1]
# for state=1:nb_states
#     if percent[state] > 0.05
#         print("Progress: ",state/nb_states*100,"%\n")
#         distances=[]
#         for step_sim=1:nb_steps
#             for carbon=1:nbC
#                 if state_matrix[step_sim,carbon] == state
#                     min_dist=V
#                     for atom=1:96
#                         dist=cell_mod.distance(traj[step_sim],cell,carbon,atom)
#                         if atom == carbon
#                             continue
#                         end
#                         if min_dist > dist
#                             min_dist = dist
#                         end
#                     end
#                     push!(distances,min_dist)
#                 end
#             end
#         end
#         min_v=0.8
#         max_v=2
#         nb_points=200
#         delta=(max_v-min_v)/nb_points
#         hist1D=zeros(nb_points)
#         for j=1:size(distances)[1]
#             for i=1:nb_points
#                 if distances[j] > min_v + (i-1)*delta && distances[j] < min_v + i*delta
#                     hist1D[i] += 1
#                 end
#             end
#         end
#         hist1D /= sum(hist1D)
#         file_distance=open(string(folder_out,"distances-",coordinances[state],"-",state,".dat"),"w")
#         for i=1:nb_points
#             write(file_distance,string(min_v+i*delta," ",hist1D[i],"\n"))
#         end
#         close(file_distance)
#     end
# end
