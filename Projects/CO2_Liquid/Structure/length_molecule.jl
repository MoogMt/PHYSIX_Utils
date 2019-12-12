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

V=8.82
T=3000

folder_base=string("/media/moogmt/Elements/CO2/",V,"/",T,"K/")

file_traj = string(folder_base,"TRAJEC_wrapped.xyz")

traj,test = filexyz.readFastFile(file_traj)
cell = cell_mod.Cell_param(V,V,V)
cut_off=1.75

nb_step=size(traj)[1]

for step=1:nb_step

    positions_local=copy(traj[step].positions)
    matrix = contact_matrix.buildMatrix( traj[step], cell, cut_off )
    molecules=graph.getGroupsFromMatrix(matrix)
    nb_molecules=size(molecules)[1]
    matrices=graph.extractAllMatrixForTrees( matrix, molecules )

    nb_atoms=size(traj[1].names)[1]
    visited=zeros(Int,nb_atoms)

    for molecule=1:nb_molecules
        if size(molecules[molecule])[1] <= 1
            continue
        end
        adjacent_molecule=getAllAdjacentVertex(matrices[molecule])
        visited=zeros(Int,size(molecules[molecule]))
        checkInfinityAndUnWrap(visited,matrices[molecule],adjacent_molecule,positions_local,cell,1,molecules[molecule],cut_off)
    end

    traj[step].positions=positions_local
end

folder_out=string(folder_base,"/Data/")
file_out=string(folder_out,"test.xyz")
filexyz.writeXYZ(file_out,traj)
