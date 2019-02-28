GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))

function buildCoordinationMatrix( traj::Vector{T1}, cell::T2, cut_off_bond::T3, max_neighbour::T4 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Real, T4 <: Int }
    # general information about the simulation
    nb_atoms=size(traj[1].names)[1]
    nb_steps=size(traj)[1]
    # Type stuff
    types=[traj[1].names[1]]
    types_number=ones(1)
    count_types=zeros(Int,1)
    for i=1:nb_atoms
        found=false
        for j=1:size(types)[1]
            if types[j] == traj[1].names[i]
                count_types[j] += 1
                found=true
            end
        end
        if ! found
            push!(types,traj[1].names[i])
            push!(types_number,size(types)[1])
            push!(count_types,1)
        end
    end
    nb_type=size(types)[1]
    type_list=zeros(Int,nb_atoms)
    for i=1:nb_atoms
        for j=1:nb_type
            if traj[1].names[i] == types[j]
                type_list[i]=j
            end
        end
    end
    # Actual computation of the coordination matrix
    coord_matrix=zeros(nb_atoms,nb_steps,max_neigh*nb_type)
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
        # Cleaning weird 3-member rings
        for atom1=1:nb_atoms-1
            for atom2=atom1+1:nb_atoms
                if bond_matrix[atom1,atom2] == 1
                    # For each bonded pair of atoms, check whether they are
                    # both bonded to another atom, forming unatural 3-member
                    # ring
                    for atom3=1:nb_atoms
                        if atom3 == atom1 || atom3 == atom1
                            continue
                        end
                        if bond_matrix[atom3,atom1] == 1 && bond_matrix[atom3,atom2] == 1
                            if cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom1,atom2) && cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom3]=0
                                bond_matrix[atom3,atom1]=0
                            elseif cell_mod.distance(traj[step_sim],cell,atom1,atom2) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom2]=0
                                bond_matrix[atom2,atom1]=0
                            else
                                bond_matrix[atom2,atom3]=0
                                bond_matrix[atom3,atom2]=0
                            end
                        end
                    end
                end
            end
        end
        # Compute coord matrix
        for atom1=1:nb_atoms
            for type=1:nb_type
                count_coord=(type-1)*max_neigh+1
                for atom2=1:nb_atoms
                    if (atom1 == atom2) || (types[type] != traj[step_sim].names[atom2])
                        continue
                    end
                    if bond_matrix[atom1,atom2] > 0
                        coord_matrix[atom1,step_sim,count_coord]=sum(bond_matrix[atom2,:])
                        count_coord += 1
                    end
                    if count_coord == type*max_neigh
                        break
                    end
                end
                #sorting by type
                for i=(type-1)*max_neigh+1:type*max_neigh-1
                    #print("type = ",type,"; i = ",i)
                    for j=i+1:type*max_neigh
                        #print(" j = ",j," ")
                        if coord_matrix[atom1,step_sim,i] < coord_matrix[atom1,step_sim,j]
                            stock = coord_matrix[atom1,step_sim,j]
                            coord_matrix[atom1,step_sim,j]=coord_matrix[atom1,step_sim,i]
                            coord_matrix[atom1,step_sim,i]=stock
                        end
                    end
                    #print("\n")
                end
            end
        end
    end
    return coord_matrix, types, type_list
end
function assignDataToStates( data::Array{T1,3}, nb_types::T2, types_number::Vector{T3} ) where { T1 <: Real, T2 <: Int , T3 <: Int }

    nb_series = size(data)[1]
    nb_steps  = size(data)[2]
    dim_data  = size(data)[3]
    state_matrix=zeros(Int, nb_series, nb_steps )

    states=zeros(Int,0,dim_data)
    record_states=zeros(Int,0)
    percent_states=zeros(Real,0)
    count_type=zeros(Real,nb_types)

    for i=1:nb_series
        print("State assignement - Progress: ",i/nb_series*100,"%\n")
        for j=1:nb_steps
            # No states, initiatilization
            if size(states)[1] == 0
                states=vcat(states,transpose(data[i,j,:]))
                count_type[types_number[i]] = count_type[types_number[i]]+ 1
                push!(record_states,types_number[i])
                push!(percent_states,1)
                state_matrix[i,j] = size(states)[1]+1
            else
                found=false
                # Loop over recorded states
                for k=1:size(states)[1]
                    # Check if types are coherent
                    if types_number[i] != record_states[k]
                        continue
                    end
                    # Compute distance between state k and data
                    dist=0
                    for l=1:dim_data
                        dist += (states[k,l]-data[i,j,l])*(states[k,l]-data[i,j,l])
                    end
                    # If dist=0 then state of data was already found
                    if dist == 0
                        count_type[types_number[i]]  = count_type[types_number[i]] + 1
                        percent_states[k] += 1
                        state_matrix[i,j] = k
                        found=true
                        break
                    end
                    # If the state was not found in database we add the state to it
                end
                if ! found
                    states=vcat(states,transpose(data[i,j,:]))
                    count_type[types_number[i]] = count_type[types_number[i]] + 1
                    push!(record_states,types_number[i])
                    push!(percent_states,1)
                    state_matrix[i,j] = size(states)[1]+1
                end
            end
        end
    end


    # Normalizing counts into percent
    for i=1:size(states)[1]
        for j=1:nb_types
            if record_states[i] == j
                percent_states[i] = percent_states[i]/count_type[j]
            end
        end
    end
    for i=1:size(states)[1]
        percent_states[i] = percent_states[i]*100
    end

    return states, percent_states, state_matrix, record_states
end
function assignDataToStates( data::Array{T1,3}, states::Array{T2,2} , nb_types::T3 , type_states::Vector{T4}, type_atoms::Vector{T5}, Err::T6 ) where { T1 <: Real, T2 <: Real , T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Bool }
    nb_series = size(data)[1]
    nb_steps  = size(data)[2]
    dim_data  = size(data)[3]
    nb_states = size(states)[1]
    state_matrix=zeros(Int, nb_series, nb_steps )
    percent_states=zeros(Real, nb_states )
    unused=0

    for j=1:nb_steps
        print("Assigning data to states - Progress: ",j/nb_steps*100,"%\n")
        for i=1:nb_series
            for l=1:nb_states
                if type_atoms[i] != type_states[l]
                    continue
                end
                d=0
                for k=1:dim_data
                    d+= ( data[i,j,k] - states[l,k])*(data[i,j,k] - states[l,k])
                end
                if d == 0
                    percent_states[l] = percent_states[l] + 1
                    state_matrix[i,j] = l
                    break
                end
            end
            if state_matrix[i,j] == -1
                if Err
                    print("Out ",i," ",j," ",dim_data," ")
                    for k=1:dim_data
                        print(data[i,j,k]," ")
                    end
                    print("\n")
                end
                unused += 1
            end
        end
    end

    # Counting the number of participant states for each type
    type_count=zeros(nb_types)
    for type=1:nb_types
        for state=1:nb_states
            if type_states[state] == type
                type_count[type] += percent_states[state]
            end
        end
    end

    # Normalizing
    for state=1:nb_states
         percent_states[state] = percent_states[state]/type_count[type_states[state]]*100
    end

    return state_matrix, percent_states
end
function isolateSignificantStates( old_states::Array{T1,2}, percent_states::Vector{T2}, cut_off_states::T3, type_states::Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real,  T4 <: Int }
    old_nb_states=size(old_states)[1]
    dim=size(old_states)[2]
    new_nb_states=0
    states_kept=zeros(0,dim)
    types_states_kept=zeros(Int,0)
    for i=1:old_nb_states
        if percent_states[i] > cut_off_states
            states_kept=[ states_kept ; transpose( old_states[i,:]) ]
            push!(types_states_kept,type_states[i])
        end
    end
    return states_kept, types_states_kept
end
function transitionMatrix( states::Array{T1,2}, state_matrix::Array{T2,2}, type_states::Vector{T3}, nb_type::T4, type_series::Vector{T5}, min_lag::T6, max_lag::T7, d_lag::T8) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real, T7 <: Int, T8 <: Int }

    nb_states=size(states)[1]
    nb_series=size(state_matrix)[1]
    nb_steps=size(state_matrix)[2]
    nb_lag_points=Int(trunc((max_lag-min_lag)/d_lag))

    count_states=zeros(Int,nb_types)
    for type=1:nb_type
        for state=1:nb_states
            if type_states[state] == type
                count_states[type] += 1
            end
        end
    end

    states_transition_matrices=[]

    for type=1:nb_type
        nb_states_loc=count_states[type]
        state_transition_probability=zeros( nb_states_loc, nb_states_loc, nb_lag_points )
        count_lag=1
        for lag=min_lag:d_lag:max_lag-1
            print("Computing Transition Matrix - Type: ",type," - Progress: ",lag/max_lag*100,"%\n")
            for serie=1:nb_series
                if type_series[serie] == type
                    for j=lag+1:nb_steps
                        if state_matrix[serie,j-lag] == 0 || state_matrix[serie,j] == 0
                            continue
                        end
                        state_transition_probability[ state_matrix[serie,j-lag]-sum(count_states[1:type-1]), state_matrix[serie,j]-sum(count_states[1:type-1]), count_lag ] += 1
                    end
                end
            end
            count_lag += 1
        end

        # Normalization
        for lag=1:nb_lag_points
            print("Normalizing Transition Matrix: ",lag/nb_lag_points*100,"%\n")
            for i=1:nb_states_loc
                sum_transition=sum( state_transition_probability[i,:,lag] )
                if sum_transition != 0
                    state_transition_probability[i,:,lag] /= sum_transition
                end
            end
        end

        if size(states_transition_matrices)[1] == 0
            states_transition_matrices=[state_transition_probability]
        else
            push!(states_transition_matrices, state_transition_probability)
        end
    end

    return states_transition_matrices
end
function chappmanKormologov( transition_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_states   = size(transition_matrix)[1]
    nb_lag_time = size(transition_matrix)[3]
    nb_lag_compare = Int(trunc(nb_lag_time/2))
    transition_matrix_kolmo= zeros(nb_states,nb_states,nb_lag_compare)
    for lag=1:nb_lag_compare
        print("Chappman Kolmogorov Test - ",lag/nb_lag_compare*100,"%\n")
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
function writeStates( file::T1 , states::Array{T2,2}, percent::Vector{T3}, types::Vector{T4}, type_list::Vector{T5} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real , T4 <: AbstractString, T5 <: Int }
    file_out=open(file,"w")
    n_dim = size( states)[2]
    nb_states=size(states)[1]
    nb_types=size(types)[1]
    write(file_out,string("  "))
    for i=1:nb_types
        for j=1:Int(n_dim/nb_types)
            write(file_out,string(types[i]," "))
        end
    end
    write(file_out,string("\n"))
    for i=1:nb_states
        write(file_out,string(types[type_list[i]]," "))
        for j=1:n_dim
            write(file_out,string(Int(states[i,j])," "))
        end
        write(file_out,string(percent[i],"\n"))
    end
    close(file_out)
    return
end
function writeStateMatrix( file::T1, state_matrix::Array{T2,2}) where { T1 <: AbstractString, T2 <: Real }
    nb_steps=size(state_matrix)[2]
    nb_series=size(state_matrix)[1]
    file_out=open(file,"w")
    for step=1:nb_steps
        write(file_out,string(step," "))
        for serie=1:nb_series
            write(file_out,string(state_matrix[serie,step]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
    return
end
function writeTransitionsMatrix( file::T1, transitions_matrix::Array{T2,3} ) where { T1 <: AbstractString, T2<: Real}
    nb_states=size(transitions_matrix)[1]
    nb_steps=size(transitions_matrix)[3]

    file_out=open(file,"w")
    for step=1:nb_steps
        write(file_out,string(step," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    return
end
function writeTransitionsMatrix( file::T1, transitions_matrix::Array{T2,3}, file2::T3, transition_matrix_CK::Array{T3,3} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    nb_states=size(transitions_matrix)[1]
    nb_steps=size(transitions_matrix)[3]

    file_out=open(file,"w")
    for step=1:nb_steps
        write(file_out,string(step," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    nb_states=size(transitions_matrix_CK)[1]
    nb_steps=size(transitions_matrix_CK)[3]

    file_out=open(file2,"w")
    for step=1:nb_steps
        write(file_out,string(step*2," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix_CK[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    return
end
function readStateMatrix( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_steps=size(lines)[1]
    nb_states=size(split(lines[1]))[1]

    state_matrix=zeros(nb_steps,nb_states)
    for step=1:nb_steps
        print("Progress: ",step/nb_steps*100,"%\n")
        line=split(lines[step])
        for states=1:nb_states
            state_matrix[step,states] = parse(Float64,line[states])
        end
    end

    return state_matrix
end
function readStates( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_states=size(lines)[1]-1
    nb_dim=size(split(lines[1]))[1]

    states=zeros(nb_states,nb_dim)
    for state=1:nb_states
        line=split(lines[state+1])
        for i=1:nb_dim
            states[state,i] = parse(Float64,line[i+1])
        end
    end

    return states, nb_states
end
function readTransitionMatrix( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_states=Int(sqrt(size(split(lines[1]) )[1]-1))
    nb_steps=size(lines)[1]

    transition_matrix=zeros(nb_states,nb_states,nb_steps)
    for step=1:nb_steps
        line=split(lines[step])
        for state_i=1:nb_states
            for state_j=1:nb_states
                transition_matrix[state_i,state_j,step] = parse(Float64,line[ (state_i-1)*nb_states+state_j+1 ])
            end
        end
    end

    return transition_matrix
end
