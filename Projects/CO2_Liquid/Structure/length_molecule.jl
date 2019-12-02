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

V=8.82
T=3000

file_traj = string(folder_in,"TRAJEC_wrapped.xyz")
folder_in = string("/media/moogmt/Elements/CO2/",V,"/",T,"K/")

traj,test = filexyz.readFastFile(file_traj)
cell = cell_mod.Cell_param(V,V,V)
cut_off=1.75

matrix = contact_matrix.buildMatrix( traj[10], cell, cut_off )
members, nb_molecules = graph.getGroupsFromMatrix( matrix )
