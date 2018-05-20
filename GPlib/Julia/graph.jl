include("contactmatrix.jl")

module graph

function searchGroupMember{ T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }( matrix::Array{T1}, list::Vector{T2}, index::T3 , group_nb::T4 )
    for i=1:size(matrix)[1]
        if matrix[index,i] > 0
            if list[i] == 0
                list[i]=group_nb
                list=searchGroupMember(matrix,list,i,group_nb)
            end
        end
    end
    return list
end

function groupsFromMatrix{ T1 <: Real, T2 <: Int }( matrix::Array{T1},  nb_vertex::T2 )
    nb_tree=0
    vertex_index=zeros(nb_vertex)
    for i=1:nb_vertex
        if vertex_index[i] == 0
            nb_tree += 1
            vertex_index=searchGroupMember(matrix,vertex_index,i,nb_tree)
        end
    end
    return vertex_index
end


    size_check=0
    size_avg=0
    for i=1:nb_mol
        write(file,string(step," ",nb_mol," "))
        size=0
        write(file,string(i," "))
        for j=1:nb_atoms
            if mol_index[j] == i
                size += 1
                write(file,string(j," "))
            end
        end
        write(file,string(size," \n"))
        push!(sizes,size)
        if size > size_check
            size_check=size
        end
        size_avg += size
    end
    push!(sizemax,size_check)
end
close(file)


end
