include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

using PyPlot

file=open("/home/moogmt/check.dat","w")

atom1=25
atom2=87

current_volume=8.82
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

file_test=open("/home/moogmt/check2.xyz","w")
write(file_test,"150\nCHECK\n")
for i=0:49
    print(string(i,"\n"))
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/",i,"_structure/")
    atoms1, cell1, ELF1 = cube_mod.readCube(string(folder,"ELF.cube"))
    pos1=cube_mod.getClosest(( atoms1.positions[atom1,:]+atoms1.positions[atom2,:])/2 ,ELF1)
    write(file,string( cell_mod.distance(atoms1,cell,atom1,atom2), " ", ELF1.matrix[ pos1[1], pos1[2], pos1[3] ] , "\n") )
end
close(file_test)

close(file)
