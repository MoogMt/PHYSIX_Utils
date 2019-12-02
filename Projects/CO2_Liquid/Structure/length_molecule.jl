GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Counts the percent of C with at least three O neighbors
# For all Temperatures and Volumes
# Needs wrapped trajectory files with thermalization
# steps removed

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
