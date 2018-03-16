include("cubefile.jl")

using PyPlot

atoms, cell, elf = cube_mod.readCube("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/0_structure/ELF.cube")

nbC=24;
nbO=[34,52,86];
