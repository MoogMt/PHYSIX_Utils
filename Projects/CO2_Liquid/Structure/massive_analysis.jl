GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"clustering.jl"))
include(string(CO2folder,"markovCO2.jl"))


# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=5001
d_lag=5
unit=0.005


V=8.82
T=3000

# for V in Volumes
#     for T in Temperatures

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
#
# if ! isfile(file)
#     continue
# end
folder_out=string(folder_in,"Data/")

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)
nb_steps=size(traj)[1]

n_neighbor=4

data=zeros(nb_steps*nbC,2*n_neighbor)
filed=open(string(folder_out,"massive_data.dat"),"w")
for step=1:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")
    for carbon=1:nbC
        distances=zeros(nbO)
        for oxygen=1:nbO
            distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
        end
        index_sort=sortperm(distances)
        data[(step-1)*nbC+carbon,1:4]=distances[index_sort[1:4]]
        distances_2nd=ones(n_neighbor)*V
        for index_=1:n_neighbor
            for carbon2=1:nbC
                if carbon2 != carbon
                    distanceOC=cell_mod.distance(traj[step],cell,carbon2,nbC+index_sort[index_])
                    if distances_2nd[index_] > distanceOC
                        distances_2nd[index_] = distanceOC
                    end
                end
            end
        end
        write(filed,string(carbon," "))
        for k=1:4
            write(filed,string(distances[index_sort[k]]," ",distances_2nd[k]," "))
            data[(step-1)*nbC+carbon,k+4] = distances_2nd[k]
        end
        write(filed,string("\n"))
    end
end
close(filed)
