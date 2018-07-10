include("xyz.jl")

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/IRENE/Pa3/333/3000K/"

atoms = filexyz.read(string(folder,"last.xyz"),1)[1]

vol=13.9875

targetVol=13.23

cell=cell_mod.Cell_param(targetVol,targetVol,targetVol)

for atom1=1:size(atoms.names)[1]
    for atom2=atom1+1:size(atoms.names)[1]
        distanceatm = cell_mod.distance(atoms,cell,atom1,atom2)
        if distanceatm < 1.0
            print("atom1: ",atom1," atom2: ",atom2," distance:",distanceatm,"\n")
        end
    end
end

Cfile=open("/home/moogmt/LargeC.cpmd","w")
Ofile=open("/home/moogmt/LargeO.cpmd","w")
for atom=1:size(atoms.names)[1]
    if atoms.names[atom] == "C"
        for i=1:3
            write(Cfile,string(atoms.positions[atom,i]," "))
        end
        write(Cfile,"\n")
    else
        for i=1:3
            write(Ofile,string(atoms.positions[atom,i]," "))
        end
        write(Ofile,"\n")
    end
end
close(Cfile)
close(Ofile)
