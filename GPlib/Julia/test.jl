include("cell.jl")
include("cubefile.jl")

using PyPlot

# Wraping atoms
atoms=cell_mod.wrap(atoms,cell)

C1=31
O1=81
O2=66
O3=70

atoms1, cell1, ELF1 = cube_mod.readCube("/home/moogmt/CO2_AIMD/ELF/0_structure/ELF.cube")
a1=cube_mod.getClosest(atoms1.positions[C1,:],ELF1)
a2=cube_mod.getClosest(atoms1.positions[O1,:],ELF1)

da=a2-a1
dir=da
dir[1]=dir[1]/dir[1]
moveMatrix=Array{Int}(7,3)
count=1
for i=0:1
    for j=0:1
        for k=0:1
            if i+j+k != 0
                moveMatrix[count,:]=[i*dir[1],j*dir[2],k*dir[3]]
                count += 1
            end
        end
    end
end

curseur=a1
list=curseur
point=curseur
dist=0

count=0
while geom.norm(curseur-a2) > 0 && count < 30
    for i=1:7
        if i==1
            point=curseur+moveMatrix[i,:]
            dist=geom.dist2Line(point,da)
        else
            dist2=geom.dist2Line(point,da)
            if dist > dist2
                point=curseur+moveMatrix[i,:]
                dist=dist2
            end
        end
    end
    print("Curseur: ", curseur[1]," ",curseur[2]," ",curseur[3]," \n")
    print("Point: ",geom.norm(curseur-a2),"\n")
    # Moving on
    curseur=point
    if count > 0
        list=[list curseur]
    end
    count++
end
