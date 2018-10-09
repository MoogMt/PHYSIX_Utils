# Loading file
include("contactmatrix.jl")

nbC=32
nbO=64

T=3000
V=9.4
cut_off=1.75

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

step=104

bonds=zeros(nbC,nbO)
for carbon=1:nbC
    for oxygen=1:nbO
        if  cell_mod.distance( traj[step], cell, carbon, oxygen+nbC ) < cut_off
            bonds[carbon,oxygen]=1
        end
    end
end

used=zeros(nbC)
for carbon1=1:nbC-1
    check=0
    for carbon2=carbon1+1:nbC
        count=0
        for oxygen=1:nbO
            if bonds[carbon1,oxygen] > 0  && bonds[carbon2,oxygen] > 0
                print("SPOTTED with C1=",carbon1," C2=",carbon2,"\n")
                count = 1
                check=1
            end
        end
        if check > 0
            break
        end
    end
end
