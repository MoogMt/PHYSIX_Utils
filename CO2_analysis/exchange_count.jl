# Loading file
include("contactmatrix.jl")

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]
Cut_Off=[1.75]

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

nbC=32
nbO=nbC*2

# for cut_off in Cut_Off
#     for T in Temperatures
#         for V in Volumes

V=9.8
T=2000
cut_off=1.75

folder_in  = string(folder_base,V,"/",T,"K/")
folder_out = string(folder_base,V,"/",T,"K/Data/")

file="TRAJEC_wrapped.xyz"
#
# if ! isfile( string(folder_in,file) )
#     continue
# end

file_out=open(string(folder_out,"exchange_count_O.dat"),"w")

traj = filexyz.readFastFile(string(folder_in,file))
cell=cell_mod.Cell_param(V,V,V)
nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]


oxygen_neighbours=zeros(nb_steps,nbO,2,2)
for step=1:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")
    for oxygen=1:nbO
        # Compute distances
        distances=zeros(nbO)
        for carbon=1:nbC
            distances[carbon]=cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
        end
        # Sorting
        for i=1:nbC
            for j=i+1:nbC
                if distances[i] > distances[j]
                    stock = distances[i]
                    distances[i] = distances[j]
                    distances[j] = stock
                end
            end
        end
        oxygen_neighbours[step,oxygen,:] = distances[1:2]
    end
end

close(file_out)

#         end
#     end
# end
