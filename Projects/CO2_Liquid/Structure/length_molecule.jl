GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz

V=8.82
T=3000

folder_in=string("/media/moogmt/Elements/CO2/",V,"/",T,"K/")

file_traj=string(folder_in,"TRAJEC_wrapped.xyz")

traj,test=readFastFile(file_traj)
