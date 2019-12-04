GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
using graph
using exp_data
using geom

V=9.1
T=3000

folder_in = string("/media/moogmt/Elements/CO2/",V,"/",T,"K/")
file_traj = string(folder_in,"TRAJEC_wrapped.xyz")

traj,test = filexyz.readFastFile(file_traj)
cell = cell_mod.Cell_param(V,V,V)
cut_off=1.75

matrix = contact_matrix.buildMatrix( traj[1], cell, cut_off )
molecules=graph.getGroupsFromMatrix(matrix)
nb_molecules=size(molecules)[1]
matrices=graph.extractAllMatrixForTrees( matrix, molecules )
adjacent_molecule=getAllAdjacentVertex(matrices[1])

start=molecules[1][1]

nb_atoms=size(molecules[1])[1]
visited=zeros(Int,nb_atoms)

function unWrapOrtho!( positions::Array{T1,2}, origin::T2, target::T3, cell::T4  ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Cell_param }
    for i=1:3
        dx=positions[origin,i]-positions[target,i]
        if abs(dx) > cell.length[i]*0.5
            signed_n = Int(trunc(dx/abs(dx)))
            n = floor(abs(dx/V))+1
            positions[target,i] += signed_n*n
        end
    end
    return
end

function recursiveExplorativeUnWrap( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}, cut_off::T8 ) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: cell_mod.Cell_param, T6 <: Int, T7 <: Int, T8 <: Real }
    visited[target]=1
    nb_neighbor=size(adjacency_table[target])[1]
    for neigh=1:nb_neighbor
        if visited[neigh] == 0
            unWrapOrtho!( positions, index_atoms[target], index_atoms[neigh], cell )
            test=recursiveExplorativeUnWrap(visited,matrix,adjacency_table,positions,cell,neigh,target,index_atoms)
            # If infinite molecule is spotted, we stop
            if ! test
                return -1
            end
        elseif neigh != target
            if geom.distance(positions[index_atoms[target],:],positions[index_atoms[neigh],:]) > cut_off
                # Spotted infinite loop; stops the search
                return -1
            end
        end
    end
    return
end

for molecule=1:nb_molecules
end


matrices=graph.extractAllMatrixForTrees(matrix,molecules)
