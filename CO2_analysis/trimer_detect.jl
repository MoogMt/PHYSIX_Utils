GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Detect and compute the number of occurences of trimers in simulations

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb

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


Volumes=[9.8]
Temperatures=[3000]

# for V in Volumes
#     for T in Temperatures

V=9.3
T=3000

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# if isfile( file )

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

nbC=32
nbO=64

cut_off = 1.8

distances_C1=zeros(nbO)
distances_C2=zeros(nbO)
index1=zeros(Int,nbO)
index2=zeros(Int,nbO)
common_neighbor=zeros(Int,2)

out_file=open(string(folder,"trimer_detect.dat"),"w")

for step=1:nb_steps
    print("V: ",V," T: ",T," Progress: ",step/nb_steps*100,"%\n")
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
                for neighbor2=neighbor1+1:5
                    if index1[neighbor1] == index2[neighbor2]
                        if cell_mod.distance(traj[step],cell,carbon1,index1[neighbor1]) < cut_off && cell_mod.distance(traj[step],cell,carbon2,index1[neighbor1]) < cut_off
                            count += 1
                            common_neighbor[count] = index1[neighbor1]
                        end
                    end
                end
            end
            if count > 1
                count_dimer += 1
                write(out_file,string(step," "))
                # print("Found at step :",step," with count of : ",count, " \n")
                write(out_file,string(carbon1," ",carbon2," ",common_neighbor[1]," ",common_neighbor[2]," ")    )
                write(out_file,string("\n"))
            end
        end
    end
    # if count_dimer > 1
    #     print("At step ",step," counted: ",count_dimer,"\n")
    # end
end
close(out_file)
#         end
#     end
# end
