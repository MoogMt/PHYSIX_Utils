include("cubefile.jl")

using PyPlot

atoms, cell, elf = cube_mod.readCube("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/0_structure/ELF.cube")

nbi=0;
nbj=0;
for i=1:32
    for j=33:size(atoms.positions)[1]
        dist=0
        for k=1:3
            dist=dist+(atoms.positions[i,k]-atoms.positions[j,k])^2
        end
        dist=sqrt(dist)
        if dist < 1.8
            nbi=i
            nbj=j
        end
        print("$i $j $dist \n")
    end
end
