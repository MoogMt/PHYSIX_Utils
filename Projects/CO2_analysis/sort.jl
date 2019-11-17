include("contactmatrix.jl")

folder="/home/moogmt/Documents/"
file=string(folder,"Super2.xyz")
traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(13.23,13.23,13.23)

nb_atoms=size(traj[1].names)[1]
atoms=traj[1]


fileC=open("/home/moogmt/C-Super.xyz","w")
fileO=open("/home/moogmt/O-Super.xyz","w")
for i=1:nb_atoms
    line=""
    for j=1:3
        line=string(line,atoms.positions[i,j]," ")
    end
    line=string(line,"\n")
    if atoms.names[i] == "C"
        write(fileC,line)
    else
        write(fileO,line)
    end
end
close(fileC)
close(fileO)
