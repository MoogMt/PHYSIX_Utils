include("cubefile.jl")

using PyPlot

# Loading file
atoms, cell, elf = cube_mod.readCube("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/0_structure/ELF.cube")

# Wraping atoms
atoms=wrap(atoms,cell)

C1=31
O1=81
O2=66
O3=70
