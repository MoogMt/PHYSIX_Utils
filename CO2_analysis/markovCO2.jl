function buildCoordinationMatrix( traj::Vector{T1}, cell::T2, cut_off_bond::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Real }
    nb_atoms=size(traj[1].names)[1]
    nb_steps=size(traj)[1]
    n_dim=9
    coord_matrix=ones(nb_steps,nbC,n_dim)*(-1)
    for step_sim=1:nb_steps

        print("Building Coordination Signal - Progress: ",step_sim/nb_steps*100,"%\n")

        # Bond Matrix
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
            for i=5:n_dim
                for j=i+1:n_dim
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
function defineStatesCoordinances()
    states=zeros(39,2)
    # CO2
    states[1,:]=[0,1]
    states[2,:]=[0,2]
    states[3,:]=[0,3]
    states[4,:]=[0,4]
    states[5,:]=[1,1]
    states[6,:]=[1,2]
    states[7,:]=[1,3]
    states[8,:]=[1,4]
    states[9,:]=[2,1]
    states[10,:]=[2,2]
    states[11,:]=[2,3]
    return states
end
function defineStatesExtendedCoordinances()
    states=zeros(86,9)
    states[1,:]=[-1,-1,-1,-1,1,-1,-1,-1,-1]
    states[2,:]=[-1,-1,-1,-1,2,-1,-1,-1,-1]
    states[3,:]=[-1,-1,-1,-1,3,-1,-1,-1,-1]
    # CO2
    states[4,:]=[-1,-1,-1,-1,1,1,-1,-1,-1]
    states[5,:]=[-1,-1,-1,-1,2,1,-1,-1,-1]
    states[6,:]=[-1,-1,-1,-1,2,2,-1,-1,-1]
    states[7,:]=[-1,-1,-1,-1,3,1,-1,-1,-1]
    states[8,:]=[-1,-1,-1,-1,3,2,-1,-1,-1]
    states[9,:]=[-1,-1,-1,-1,3,3,-1,-1,-1]
    # CO3
    states[10,:]=[-1,-1,-1,-1,2,2,2,-1,-1]
    states[11,:]=[-1,-1,-1,-1,2,2,1,-1,-1]
    states[12,:]=[-1,-1,-1,-1,2,1,1,-1,-1]
    states[13,:]=[-1,-1,-1,-1,1,1,1,-1,-1]
    #
    states[14,:]=[-1,-1,-1,-1,3,1,1,-1,-1]
    states[15,:]=[-1,-1,-1,-1,3,2,1,-1,-1]
    states[16,:]=[-1,-1,-1,-1,3,2,2,-1,-1]
    states[17,:]=[-1,-1,-1,-1,3,3,1,-1,-1]
    states[18,:]=[-1,-1,-1,-1,3,3,2,-1,-1]
    states[19,:]=[-1,-1,-1,-1,3,3,3,-1,-1]
    # CO4
    states[20,:]=[-1,-1,-1,-1,2,2,2,2,-1]
    states[21,:]=[-1,-1,-1,-1,2,2,2,1,-1]
    states[22,:]=[-1,-1,-1,-1,2,2,1,1,-1]
    states[23,:]=[-1,-1,-1,-1,2,1,1,1,-1]
    states[24,:]=[-1,-1,-1,-1,1,1,1,1,-1]
    states[25,:]=[-1,-1,-1,-1,3,1,1,1,-1]
    states[26,:]=[-1,-1,-1,-1,3,2,1,1,-1]
    states[27,:]=[-1,-1,-1,-1,3,2,2,1,-1]
    states[28,:]=[-1,-1,-1,-1,3,2,2,2,-1]
    states[29,:]=[-1,-1,-1,-1,3,3,1,1,-1]
    states[30,:]=[-1,-1,-1,-1,3,3,2,1,-1]
    states[31,:]=[-1,-1,-1,-1,3,3,2,2,-1]
    states[32,:]=[-1,-1,-1,-1,3,3,3,1,-1]
    states[33,:]=[-1,-1,-1,-1,3,3,3,2,-1]
    states[34,:]=[-1,-1,-1,-1,3,3,3,3,-1]
    # CO5
    states[35,:]=[-1,-1,-1,-1,2,2,2,2,2]
    states[36,:]=[-1,-1,-1,-1,2,2,2,2,1]
    states[37,:]=[-1,-1,-1,-1,2,2,2,1,1]
    states[38,:]=[-1,-1,-1,-1,2,2,1,1,1]
    states[39,:]=[-1,-1,-1,-1,2,1,1,1,1]
    # CCO1
    states[40,:]=[2,-1,-1,-1,2,-1,-1,-1,-1]
    states[41,:]=[2,-1,-1,-1,1,-1,-1,-1,-1]
    states[42,:]=[3,-1,-1,-1,2,-1,-1,-1,-1]
    states[43,:]=[3,-1,-1,-1,1,-1,-1,-1,-1]
    states[44,:]=[4,-1,-1,-1,2,-1,-1,-1,-1]
    states[45,:]=[4,-1,-1,-1,1,-1,-1,-1,-1]
    # CCO2
    states[46,:]=[2,-1,-1,-1,2,2,-1,-1,-1]
    states[47,:]=[2,-1,-1,-1,2,1,-1,-1,-1]
    states[48,:]=[2,-1,-1,-1,1,1,-1,-1,-1]
    states[49,:]=[3,-1,-1,-1,2,2,-1,-1,-1]
    states[50,:]=[3,-1,-1,-1,2,1,-1,-1,-1]
    states[51,:]=[3,-1,-1,-1,1,1,-1,-1,-1]
    states[52,:]=[4,-1,-1,-1,2,2,-1,-1,-1]
    states[53,:]=[4,-1,-1,-1,2,1,-1,-1,-1]
    states[54,:]=[4,-1,-1,-1,1,1,-1,-1,-1]
    # CCO3
    states[55,:]=[2,-1,-1,-1,2,2,2,-1,-1]
    states[56,:]=[2,-1,-1,-1,2,2,1,-1,-1]
    states[57,:]=[2,-1,-1,-1,2,1,1,-1,-1]
    states[58,:]=[2,-1,-1,-1,1,1,1,-1,-1]
    states[59,:]=[3,-1,-1,-1,2,2,2,-1,-1]
    states[60,:]=[3,-1,-1,-1,2,2,1,-1,-1]
    states[61,:]=[3,-1,-1,-1,2,1,1,-1,-1]
    states[62,:]=[3,-1,-1,-1,1,1,1,-1,-1]
    states[63,:]=[4,-1,-1,-1,2,2,2,-1,-1]
    states[64,:]=[4,-1,-1,-1,2,2,1,-1,-1]
    states[65,:]=[4,-1,-1,-1,2,1,1,-1,-1]
    states[66,:]=[4,-1,-1,-1,1,1,1,-1,-1]
    #
    states[67,:]=[2,-1,-1,-1,2,2,2,1,-1]
    states[68,:]=[2,-1,-1,-1,2,2,1,1,-1]
    states[69,:]=[2,-1,-1,-1,2,1,1,1,-1]
    states[70,:]=[2,-1,-1,-1,1,1,1,1,-1]
    states[71,:]=[3,-1,-1,-1,2,2,2,1,-1]
    states[72,:]=[3,-1,-1,-1,2,2,1,1,-1]
    states[73,:]=[3,-1,-1,-1,2,1,1,1,-1]
    states[74,:]=[3,-1,-1,-1,1,1,1,1,-1]
    states[75,:]=[4,-1,-1,-1,2,2,2,1,-1]
    states[76,:]=[4,-1,-1,-1,2,2,1,1,-1]
    states[77,:]=[4,-1,-1,-1,2,1,1,1,-1]
    states[78,:]=[4,-1,-1,-1,1,1,1,1,-1]
    #
    states[79,:]=[2,-1,-1,-1,2,2,2,2,-1]
    states[80,:]=[2,-1,-1,-1,2,2,2,1,-1]
    states[81,:]=[2,-1,-1,-1,2,2,1,1,-1]
    states[82,:]=[2,-1,-1,-1,2,1,1,1,-1]
    states[83,:]=[3,-1,-1,-1,2,2,2,2,-1]
    states[84,:]=[3,-1,-1,-1,2,2,2,1,-1]
    states[85,:]=[3,-1,-1,-1,2,2,1,1,-1]
    states[86,:]=[3,-1,-1,-1,2,1,1,1,-1]
    return states
end
function defineStatesExtendedCoordinancesO()
    states=zeros(39,8)
    # CO2
    states[1,:]=[2,-1,-1]
    states[1,:]=[3,-1,-1]
    states[1,:]=[4,-1,-1]
    states[1,:]=[2,2,-1]
    states[1,:]=[2,3,-1]
    states[1,:]=[2,4,-1]
    states[1,:]=[3,2,-1]
    states[1,:]=[3,3,-1]
    states[1,:]=[3,4,-1]
    states[1,:]=[4,2,-1]
    states[1,:]=[4,3,-1]
    states[1,:]=[4,4,-1]
    #
    states[1,:]=[2,2,2]
    states[1,:]=[2,3,2]
    states[1,:]=[2,4,2]
    states[1,:]=[3,2,2]
    states[1,:]=[3,3,2]
    states[1,:]=[3,4,2]
    states[1,:]=[4,2,2]
    states[1,:]=[4,3,2]
    states[1,:]=[4,4,2]
    #
    states[1,:]=[2,2,3]
    states[1,:]=[2,3,3]
    states[1,:]=[2,4,3]
    states[1,:]=[3,2,3]
    states[1,:]=[3,3,3]
    states[1,:]=[3,4,3]
    states[1,:]=[4,2,3]
    states[1,:]=[4,3,3]
    states[1,:]=[4,4,3]
    #
    states[1,:]=[2,2,4]
    states[1,:]=[2,3,4]
    states[1,:]=[2,4,4]
    states[1,:]=[3,2,4]
    states[1,:]=[3,3,4]
    states[1,:]=[3,4,4]
    states[1,:]=[4,2,4]
    states[1,:]=[4,3,4]
    states[1,:]=[4,4,4]
    return states
end
function assignDataToStates( data::Array{T1,3}, states::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_data_point=size(data)[1]
    nb_series = size(data)[2]
    dim_data = size(data)[3]
    nb_states = size(states)[1]
    state_matrix=ones(Int, nb_data_point, nb_series )*(-1)
    count_states=zeros( nb_states )
    unused=0
    for i=1:nb_data_point
        print("Assigning data to states - Progress: ",i/nb_data_point*100,"%\n")
        for j=1:nb_series
            for l=1:nb_states
                d=0
                for k=1:dim_data
                    d+= ( data[i,j,k] - states[l,k])*(data[i,j,k] - states[l,k])
                end
                if d == 0
                    count_states[l] += 1
                    state_matrix[i,j] = l
                    break
                end
            end
            if state_matrix[i,j] == -1
                print("test ")
                for k=1:dim_data
                    print(data[i,j,k]," ")
                end
                print("\n")
                unused += 1
            end
        end
    end
    return state_matrix, count_states/sum(count_states)*100, unused/(unused+sum(count_states))*100
end
function isolateSignificantStates( old_states::Array{T1,2}, percent_states::Vector{T2}, cut_off_states::T3 ) where { T1 <: Real, T2 <: Real , T3 <: Real }
    old_nb_states=size(old_states)[1]
    dim=size(old_states)[2]
    new_nb_states=0
    states_kept=zeros(0,dim)
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
                if  state_matrix[j-lag,i] == -1 ||  state_matrix[j,i] == -1
                    continue
                end
                states_transition_probability[ state_matrix[j-lag,i], state_matrix[j,i], count_lag ] += 1
            end
        end
        count_lag += 1
    end

    # Normalization
    for lag=1:nb_lag_points
        for i=1:nb_states
            states_transition_probability[i,:,lag] /= sum( states_transition_probability[i,:,lag] )
        end
    end

    return states_transition_probability
end
function chappmanKormologov( transition_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_states   = size(transition_matrix)[1]
    nb_lag_time = size(transition_matrix)[3]
    nb_lag_compare = Int(trunc(nb_lag_time/2))
    transition_matrix_kolmo= zeros(nb_states,nb_states,nb_lag_compare)
    for lag=1:nb_lag_compare
        for i=1:nb_states
            for j=1:nb_states
                for k=1:nb_states
                    transition_matrix_kolmo[i,j,lag] += transition_matrix[i,k,lag]*transition_matrix[k,j,lag]
                end
            end
        end
    end
    return transition_matrix_kolmo
end
function writeStates( file::T1 , states::Array{T2,2}, percent::Vector{T3}) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    file_out=open(file,"w")
    n_dim = size( states)[2]
    nb_states=size(states)[1]
    for i=1:nb_states
        for j=1:n_dim
            write(file_out,string(states[i,j]," "))
        end
        write(file_out,string(percent[i],"\n"))
    end
    close(file_out)
    return
end
