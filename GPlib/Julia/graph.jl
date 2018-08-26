module graph_mod

export searchGroupMember, groupsFromMatrix, getSizeTree

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
            vertex_index = searchGroupMember(matrix,vertex_index,i,nb_tree)
        end
    end
    return nb_tree, vertex_index
end

function getSizeTrees{ T1 <: Int }( vertex_index::Vector{T1} )
    sizes=[]
    max=0
    tree_index=unique(vertex_index)
    for tree in tree_index
        size=0
        for i=1:size(vertex_index)[1]
            if tree == vertex_index
                size+=1
            end
        end
        if size > max
            max=size
        end
        push!(sizes,size)
    end
    return sizes, max
end

end
