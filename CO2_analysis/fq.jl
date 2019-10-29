GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

T=3000
V=8.82

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

nbO=64
nbC=32

nb_nn=8

data_file=string(folder_out,"data_CO-",nb_nn,".dat")

dataCO=[]

erase=true

if ! isfile(data_file) || erase

    traj=filexyz.readFastFile(file)
    cell=cell_mod.Cell_param(V,V,V)

    nb_steps=size(traj)[1]

    dataCO=zeros(nb_steps,nbC,nb_nn)
    for t=1:nb_steps
        for carbon=1:nbC
            distances=zeros(nbO)
            for oxygen=1:nbO
                distances[oxygen] = cell_mod.distance(traj[t],cell,carbon,nbC+oxygen)
            end
            sort!(distances)
            dataCO[t,carbon,:]=distances[1:nb_nn]
        end
    end

    file_out=open(data_file,"w")
    for i=1:nb_steps
        for j=1:nbC
            for k=1:nb_nn
                write(file_out,string(dataCO[i,j,k]," "))
            end
            write(file_out,string("\n"))
        end
    end
    close(file_out)
else

    file_in=open(data_file)
    lines=readlines(file_in)
    close(file_in)

    nb_steps=Int(size(lines)[1]/nbC)

    dataCO=zeros(nb_steps,nbC,nb_nn)
    for i=1:nb_steps
        for j=1:nbC
            for k=1:nb_nn
                dataCO[i,j,k] = parse(Float64,split(lines[(i*nbC-1)+j])[k])
            end
        end
    end
