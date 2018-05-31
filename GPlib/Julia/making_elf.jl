include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

using PyPlot

C1=31
O1=81
O2=66
O3=70

file=open("/home/moogmt/check.dat","w")

current_volume=8.82
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/",0,"_structure/")
atoms1, cell1, ELF1 = cube_mod.readCube(string(folder,"ELF.cube"))

for i=1:96
    for j=i+1:96
        if cell_mod.distance( atoms1, cell, i, j ) < 1.8
            print(string(i," ", j, "\n"))
        end
    end
end

for i=0:49
    print(string(i,"\n"))
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/",i,"_structure/")
    atoms1, cell1, ELF1 = cube_mod.readCube(string(folder,"ELF.cube"))
    pos1=cube_mod.getClosest((atoms1.positions[1,:]+atoms1.positions[59,:])/2-ELF1.origin,ELF1)
    write(file,string( cell_mod.distance(atoms1,cell,1,59), " ", ELF1.matrix[ pos1[1], pos1[2], pos1[3] ] , "\n") )
end

close(file)
