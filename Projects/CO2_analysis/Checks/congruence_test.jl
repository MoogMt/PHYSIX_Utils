GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using cpmd

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# T,V
V=10.0
T=3000

# Input folder
folder_in=string(folder_base,V,"/",T,"K/1-run/")
file_input=string(folder_in,"input")

time_step=cpmd.readInputTimestep( file_input )
stride_stress=cpmd.readIntputStrideStress( file_input )
stride_traj=cpmd.readIntputStrideTraj(file_input)
