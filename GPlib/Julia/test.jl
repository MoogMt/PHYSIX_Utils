include("cell.jl")
include("cubefile.jl")

using PyPlot

# Wraping atoms
atoms=cell_mod.wrap(atoms,cell)

C1=31
O1=81
O2=66
O3=70

atoms, cell, ELF = cube_mod.readCube("/home/moogmt/CO2_AIMD/ELF/0_structure/ELF.cube")
