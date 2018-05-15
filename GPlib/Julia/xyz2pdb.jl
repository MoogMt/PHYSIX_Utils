#Adapt to target the lib_directory
folder_script="/home/moogmt/PHYSIX_Utils/GPlib/Julia/"

# Including necessary files...
include("xyz.jl")
include("pdb.jl")

# Getting path of file
file_xyz="/home/moogmt/9.8/position.xyz"

# path to write
file_out="/home/moogmt/large.pdb"

# Creating cell
cell=cell_mod.Cell_param(9.325,9.325,9.325)

# Reading XYZ file
atoms = filexyz.readFastFile(file_xyz)[1]

pdb.writeStep(atoms,cell,file_out)
