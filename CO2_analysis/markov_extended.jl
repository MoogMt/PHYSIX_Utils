GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))

function buildCoordinationMatrix( traj::Vector{T1}, cell::T2, cut_off_bond::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Real }
    nb_atoms=size(traj[1].names)[1]
    nb_steps=size(traj)[1]
    coord_matrix=ones(nb_steps,nbC,8)*(-1)
    for step_sim=1:nb_steps
        print("Building Coordination Signal - Progress: ",step_sim/nb_steps*100,"%\n")
        bond_matrix=zeros(nb_atoms,nb_atoms)
        for atom1=1:nb_atoms
            for atom2=atom1+1:nb_atoms
                if cell_mod.distance( traj[step_sim], cell, atom1, atom2) < cut_off_bond
                    bond_matrix[atom1,atom2]=1
                    bond_matrix[atom2,atom1]=1
                end
            end
        end

        for carbon=1:nbC
            count_coord=1
            for carbon2=1:nbC
                if carbon == carbon2
                    continue
                end
                if bond_matrix[carbon,carbon2] > 0
                    coord_matrix[step_sim,carbon,count_coord]=sum(bond_matrix[carbon2,:])
                    count_coord += 1
                end
            end
            count_coord=1
            for oxygen=1:nbO
                if bond_matrix[carbon,nbC+oxygen] > 0
                    if count_coord > 4
                        continue
                    end
                    coord_matrix[step_sim,carbon,4+count_coord]=sum(bond_matrix[nbC+oxygen,:])
                    count_coord += 1
                end
            end
            # sort
            for i=1:4
                for j=i+1:4
                    if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                        stock=coord_matrix[step_sim,carbon,i]
                        coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                        coord_matrix[step_sim,carbon,j]=stock
                    end
                end
            end
            for i=5
                for j=i+1:8
                    if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                        stock=coord_matrix[step_sim,carbon,i]
                        coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                        coord_matrix[step_sim,carbon,j]=stock
                    end
                end
            end
        end
    end
    return coord_matrix
end

function defineStates()
    states=zeros(39,8)
    # CO2
    states[1,:]=[-1,-1,-1,-1,1,1,-1,-1]
    states[2,:]=[-1,-1,-1,-1,2,1,-1,-1]
    states[3,:]=[-1,-1,-1,-1,2,2,-1,-1]
    # CO3
    states[4,:]=[-1,-1,-1,-1,2,2,2,-1]
    states[5,:]=[-1,-1,-1,-1,2,2,1,-1]
    states[6,:]=[-1,-1,-1,-1,2,1,1,-1]
    states[7,:]=[-1,-1,-1,-1,1,1,1,-1]
    # CO4
    states[8,:]=[-1,-1,-1,-1,2,2,2,2]
    states[9,:]=[-1,-1,-1,-1,2,2,2,1]
    states[10,:]=[-1,-1,-1,-1,2,2,1,1]
    states[11,:]=[-1,-1,-1,-1,2,1,1,1]
    states[12,:]=[-1,-1,-1,-1,1,1,1,1]
    # CCO1
    states[13,:]=[2,-1,-1,-1,2,-1,-1,-1]
    states[14,:]=[2,-1,-1,-1,1,-1,-1,-1]
    states[15,:]=[3,-1,-1,-1,2,-1,-1,-1]
    states[16,:]=[3,-1,-1,-1,1,-1,-1,-1]
    states[17,:]=[4,-1,-1,-1,2,-1,-1,-1]
    states[18,:]=[4,-1,-1,-1,1,-1,-1,-1]
    # CCO2
    states[19,:]=[2,-1,-1,-1,2,2,-1,-1]
    states[20,:]=[2,-1,-1,-1,2,1,-1,-1]
    states[21,:]=[2,-1,-1,-1,1,1,-1,-1]
    states[22,:]=[3,-1,-1,-1,2,2,-1,-1]
    states[23,:]=[3,-1,-1,-1,2,1,-1,-1]
    states[24,:]=[3,-1,-1,-1,1,1,-1,-1]
    states[25,:]=[4,-1,-1,-1,2,2,-1,-1]
    states[26,:]=[4,-1,-1,-1,2,1,-1,-1]
    states[27,:]=[4,-1,-1,-1,1,1,-1,-1]
    # CCO3
    states[28,:]=[2,-1,-1,-1,2,2,2,-1]
    states[29,:]=[2,-1,-1,-1,2,2,1,-1]
    states[30,:]=[2,-1,-1,-1,2,1,1,-1]
    states[31,:]=[2,-1,-1,-1,1,1,1,-1]
    states[32,:]=[3,-1,-1,-1,2,2,2,-1]
    states[33,:]=[3,-1,-1,-1,2,2,1,-1]
    states[34,:]=[3,-1,-1,-1,2,1,1,-1]
    states[35,:]=[3,-1,-1,-1,1,1,1,-1]
    states[36,:]=[4,-1,-1,-1,2,2,2,-1]
    states[37,:]=[4,-1,-1,-1,2,2,1,-1]
    states[38,:]=[4,-1,-1,-1,2,1,1,-1]
    states[39,:]=[4,-1,-1,-1,1,1,1,-1]
    return states
end

function assignDataToStates( data::Array{T1,3}, states::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_data_point=size(data)[1]
    nb_series = size(data)[2]
    dim_data = size(data)[3]
    nb_states = size(states)[1]
    state_matrix=zeros(Int, nb_data_point, nb_series )
    count_states=zeros( nb_states )
    for i=1:nb_data_point
        print("Assigning data to states - Progress: ",i/nb_data_point*100,"%\n")
        for j=1:nb_series
            index=1
            d=0
            for k=1:dim_data
                d+= ( data[i,j,k] - states[1,k])*(data[i,j,k] - states[1,k])
            end
            d_min=d
            for k=2:nb_states
                d=0
                for l=1:dim_data
                    d+= ( data[i,j,l] - states[k,l])*(data[i,j,l] - states[k,l])
                end
                if d < d_min
                    index=k
                    d_min=d
                end
            end
            count_states[index] += 1
            state_matrix[i,j]=index
        end
    end
    return state_matrix, count_states/sum(count_states)*100
end

function isolateSignificantStates( old_states::Array{T1,2}, percent_states::Vector{T2}, cut_off_states::T3 ) where { T1 <: Real, T2 <: Real , T3 <: Real }
    old_nb_states=size(old_states)[1]
    new_nb_states=0
    states_kept=zeros(0,8)
    for i=1:old_nb_states
        if percent_states[i] > cut_off_states
            states_kept=[ states_kept ; transpose( old_states[i,:]) ]
        end
    end
    return states_kept
end

function transitionMatrix( states::Array{T1,2}, state_matrix::Array{T2,2}, min_lag::T3, max_lag::T4, d_lag::T5) where { T1 <: Real, T2 <: Real, T3 <: Real, T4<:Int, T5 <: Int }

    nb_states=size(states)[1]
    nb_data_point=size(state_matrix)[1]
    nb_series = size(state_matrix)[2]
    nb_lag_points=Int(trunc((max_lag-min_lag)/d_lag))

    states_transition_probability=zeros(Float64,nb_states,nb_states,nb_lag_points)

    # Chappman Kolmogorov test
    count_lag=1
    for lag=min_lag:d_lag:max_lag-1
        print("Chappman Kolmogorov Test - Progress: ",lag/max_lag*100,"%\n")
        for i=1:nb_series
            for j=lag+1:nb_data_point
                states_transition_probability[ state_matrix[j-lag,i], state_matrix[j,i], count_lag ] += 1
            end
        end
        count_lag += 1
    end

    # Normalization
    for lag=1:nb_lag_points
        for i=1:nb_states
            states_transition_probability[:,i,lag] /= sum( states_transition_probability[:,i,lag] )
        end
    end

    return states_transition_probability
end

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1
min_lag=1
max_lag=2001
d_lag=5

V=8.82
T=3000

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

states=defineStates()
data=buildCoordinationMatrix( filexyz.readFastFile(file),  cell_mod.Cell_param(V,V,V), cut_off_bond )
case_matrix, percent = assignDataToStates( data , states )

statistic_states = isolateSignificantStates( states, percent, cut_off_states )
state_matrix, percent = assignDataToStates( data , statistic_states )
transition_matrix = transitionMatrix( statistic_states, state_matrix, min_lag, max_lag, d_lag )

#
# cases_keep=[]
# count_cases_keep=[]
# case_matrix_keep=[]
# case_proba=[]
# for V in Volumes
#     for T in Temperatures
#         for cut_off in Cut_Off
#
#
#             #
#         end
#     end
# end
#
# file_cases_keep=open(string(folder_out,"cases-",cut_off,".dat"),"w")
# if size(count_cases_keep)[1] == 1
#     for i=1:8
#         write(file_cases_keep,string(cases_keep[i]," "))
#     end
#     write(file_cases_keep,string(count_cases_keep,"\n"))
# else
#     for i=size(count_cases_keep)[1]
#         for j=1:8
#             write(file_cases_keep,string(cases_keep[i,j]," "))
#         end
#         write(file_cases_keep,string(count_cases_keep,"\n"))
#     end
# end
# close(file_cases_keep)
#
#
# file_proba=open(string(folder_out,"proba_transition-",cut_off,".dat"),"w")
# for lag=1:size(case_proba)[3]
#     if size(count_cases_keep)[1] == 1
#         write(file_proba,string(case_proba[1,1,lag]," "))
#         continue
#     end
#     for i=1:size(case_proba)[1]
#         for j=1:size(case_proba)[2]
#             write(file_proba,string(case_proba[i,j,lag]," "))
#         end
#     end
#     write(file_proba,string("\n"))
# end
# close(file_proba)
#
# file_lag=open(string(folder_out,"markovian_evolution_test-",cut_off,".dat"),"w")
# file_out_2=open(string(folder_out,"markovian_evolution_test2-",cut_off,".dat"),"w")
# for lag=1:count_lag-1
#     d_test=0
#     d_test2=0
#     count2=0
#     for i=1:size(cases_keep)[1]
#         for j=1:size(cases_keep)[1]
#             # Specific i-j transition
#             p_truth=case_proba[i,j,lag]
#             p_test=0
#             p_test2=0
#             for k=1:size(cases_keep)[1]
#                 p_test += case_proba[i,k,Int(trunc(lag/2))+1]*case_proba[k,j,Int(trunc(lag/2))+1]
#             end
#             d_test += abs(p_truth-p_test)
#             d_test2 += d_test*d_test
#             count2 += 1
#             write(file_lag,string(p_truth," ",p_test," "))
#         end
#     end
#     write(file_lag,string("\n"))
#     write(file_out_2,string(lag*d_lag*0.005," ",d_test/count2," ",d_test2/count2-(d_test/count2)^2,"\n"))
# end
# close(file_lag)
# close(file_out_2)
#
# file_lag=open(string(folder_out,"markovian_evolution_test-3_",cut_off,".dat"),"w")
# for lag=2:count_lag-1
#     case_proba_test=case_proba[:,:,1]^lag
#     write(file_lag,string(lag*d_lag*0.005," "))
#     for i=1:size(cases_keep)[1]
#         for j=1:size(cases_keep)[1]
#             # Specific i-j transition
#             write(file_lag,string(case_proba[i,j,lag]," ",case_proba_test[i,j]," "))
#         end
#     end
#     write(file_lag,string("\n"))
# end
#
# file_cases=open(string(folder_out,"cases_occurences-",cut_off,".dat"),"w")
# for step_sim=1:nb_steps
#     for carbon=1:nbC
#         write(file_cases,string(case_matrix_keep[step_sim,carbon]," "))
#     end
#     write(file_cases,string("\n"))
# end
# close(file_cases)
