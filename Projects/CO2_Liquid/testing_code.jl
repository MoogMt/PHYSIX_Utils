GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using cpmd

folder_in=string("/media/moogmt/Elements/CO2/8.82/3000K/1-run/")

file=string(folder_in,"FTRAJECTORY")

positions,velocity,forces=cpmd.readFTRAJ(file)
