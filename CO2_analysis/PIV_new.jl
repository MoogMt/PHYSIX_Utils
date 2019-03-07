GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

T=3000
V=8.82

n_run=1

folder_in=string(folder_base,V,"/",T,"K/",n_run,"-run/")

print("Computing Data\n")
traj=filexyz.readFastFile(string(folder_in,"TRAJEC.xyz"))
cell=cell_mod.Cell_param(V,V,V)

traj=traj[1:100]
nb_structure=size(traj)[1]
nb_atoms=size(traj[1].names)[1]
nbC=32
nbO=64

# Wrapping
for i=1:nb_structure
    cell_mod.wrap(traj[i],cell)
end

# ComputePIV
d0=1.8
m=10
n=4
size_piv=Int(nbC*nbO)
piv=zeros(size_piv,nb_structure)
for step=1:nb_structure
    print("PIV computation progress: ",step/nb_structure*100,"%\n")
    count=1
    for carbon=1:nbC
        for oxygen=1:nbO
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,nbC+oxygen),d0,n,m)
            count = count + 1
        end
    end
    # Sort
    for atom1=1:size_piv-1
        for atom2=atom1+1:size_piv
            if piv[atom1,step] < piv[atom2,step]
                stock=piv[atom1,step]
                piv[atom1,step] = piv[atom2,step]
                piv[atom2,step] = stock
            end
        end
    end
end

nb_box=1000
box_piv=zeros(nb_box,nb_structure)
for structure=1:nb_structure
    print("PIV computation progress: ",structure/nb_structure*100,"%\n")
    for carbon=1:nbC
        for oxygen=1:nbO
            dist=utils.switchingFunction(cell_mod.distance(traj[structure],cell,carbon,nbC+oxygen),d0,n,m)
            for box=1:nb_box
                if (box-1)/nb_box < dist && box/nb_box > dist
                    box_piv[1:box,structure] = box_piv[1:box,structure].+1
                    break
                end
            end
        end
    end
end

file_out=open(string("/home/moogmt/hist_piv-",nb_box,".dat"),"w")
for i=1:nb_box
    write(file_out,string(i/nb_box," "))
    for structure=1:nb_structure
        write(file_out,string(box_piv[i,structure]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

file_out=open(string("/home/moogmt/piv-",nb_box,".dat"),"w")
for i=1:size_piv
    write(file_out,string(i," "))
    for structure=1:nb_structure
        write(file_out,string(piv[i,structure]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

file_out=open(string("/home/moogmt/test-",nb_box,".dat"),"w")
for structure1=1:nb_structure
    for structure2=1:nb_structure
        dist_piv=0
        for i=1:size_piv
            dist_piv += abs(piv[i,structure1]-piv[i,structure2])
        end
        dist_hist=0
        for i=1:nb_box
            dist_hist += abs( (box_piv[i,structure1]-box_piv[i,structure2])/nb_box )
        end
        write(file_out,string(dist_piv," ",dist_hist,"\n"))
    end
end
close(file_out)
