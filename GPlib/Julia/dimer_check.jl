# Loading file
include("contactmatrix.jl")

function sortIndex( indexes, x )
    sizex=size(x)[1]
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stock2=indexes[i]
                indexes[i]=indexes[j]
                x[i]=x[j]
                indexes[j]=stock2
                x[j]=stock
            end
        end
    end
    return
end


nbC=32
nbO=64

distances_C1=zeros(nbO)
distances_C2=zeros(nbO)
index1=zeros(Int,nbO)
index2=zeros(Int,nbO)

stop_100=20000

common_neighbor=zeros(Int,2)

V=9.4
T=3000
cut_off=1.75

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]


for step=1:stop_100
    if step > nb_steps
        break
    end

    count_dimer = 0
    for carbon1 = 1:nbC

        for oxygen = 1:nbO
            distances_C1[oxygen] = cell_mod.distance(traj[step],cell,carbon1,oxygen)
        end
        for i=1:nbO
            index1[i]=i
        end
        sortIndex(index1, distances_C1)
        for carbon2 = carbon1+1:nbC
            for oxygen = 1:nbO
                distances_C2[oxygen] = cell_mod.distance(traj[step],cell,carbon2,oxygen)
            end
            for i=1:nbO
                index2[i]=i
            end
            sortIndex(index2, distances_C2)
            count = 0
            for neighbor1=2:5
                for neighbor2=2:5
                    if index1[neighbor1] == index2[neighbor2]
                        if cell_mod.distance(traj[step],cell,carbon1,index1[neighbor1]) < cut_off
                            if  cell_mod.distance(traj[step],cell,carbon2,index2[neighbor2]) < cut_off
                                count += 1
                                common_neighbor[count] = index1[neighbor1]
                            end
                        end
                    end
                end
            end
            if count > 1
                count_dimer += 1
                print(out_file,string(step," "))
                print(out_file,string(carbon1," ",carbon2," ",common_neighbor[1]," ",common_neighbor[2]," ")    )
                print(out_file,string("\n"))
            end
        end
    end
end
