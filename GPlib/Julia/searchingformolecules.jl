# Loading file
include("contactmatrix.jl")

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
file=string(folder,"TRAJEC_wrapped.xyz")

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(8.82,8.82,8.82)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Stride
stride=5
unit=0.005

for step=1:nb_steps
    
end
