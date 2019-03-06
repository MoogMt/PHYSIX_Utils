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
piv=zeros(Int(nb_atoms*(nb_atoms-1)/2),nb_structure)
for step=1:nb_structure
    print("PIV computation progress: ",step/nb_structure*100,"%\n")
    count=1
    start=count
    for carbon=1:nbC
        for carbon2=carbon+1:nbC
            piv[count,step] = cell_mod.distance(traj[step],cell,carbon,carbon2)
            count = count + 1
        end
    end
    # # Sort
    for atom=start:count-2
        for atom2=atom+1:count-1
            if piv[atom,step] > piv[atom2,step]
                stock=piv[atom,step]
                piv[atom,step] = piv[atom2,step]
                piv[atom2,step] = stock
            end
        end
    end
    start=count
    for oxygen=1:nbO
        for oxygen2=oxygen+1:nbO
            piv[count,step] = cell_mod.distance(traj[step],cell,nbC+oxygen,nbC+oxygen2)
            count = count + 1
        end
    end
    # Sort
    for atom=start:count-2
        for atom2=atom+1:count-1
            if piv[atom,step] > piv[atom2,step]
                stock=piv[atom,step]
                piv[atom,step] = piv[atom2,step]
                piv[atom2,step] = stock
            end
        end
    end
    start=count
    for carbon=1:nbC
        for oxygen=1:nbO
            piv[count,step] = cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
            count = count + 1
        end
    end
    # Sort
    for atom=start:count-2
        for atom2=atom+1:count-1
            if piv[atom,step] > piv[atom2,step]
                stock=piv[atom,step]
                piv[atom,step] = piv[atom2,step]
                piv[atom2,step] = stock
            end
        end
    end
end


d0=1.8
n=4
m=10

file_out=open(string("/home/moogmt/PIV-",d0,"-",n,"-",m,".dat"),"w")
for i=1:size(piv)[1]
    print("Progress : ",i/size(piv)[1]*100,"%\n")
    write(file_out,string(i," "))
    for j=1:size(piv)[2]
        write(file_out,string(utils.switchingFunction(piv[i,j],d0,n,m)," "))
    end
    write(file_out,string("\n"))
end
close(file_out)
distances=zeros(nb_structure,nb_structure)
for step1=1:nb_structure
    print("Computing distance matrix: ",step1/nb_structure*100,"%\n")
    for step2=step1+1:nb_structure
        for i=1:size(piv)[1]
            distances[step1,step2] = distances[step1,step2] + (utils.switchingFunction(piv[i,step1],d0,n,m)-utils.switchingFunction(piv[i,step2],d0,n,m))*(utils.switchingFunction(piv[i,step1],d0,n,m)-utils.switchingFunction(piv[i,step2],d0,n,m))
        end
        distances[step1,step2] = sqrt(distances[step1,step2])
        distances[step2,step1] = distances[step1,step2]
    end
end

nb_box=200
hist_PIV=zeros(nb_box,nb_structure)
for structure=1:nb_structure
    for
    for box=1:nb_box

    end
end


file_energy=open(string(folder_in,"EKS_base"))
lines=readlines(file_energy)
close(file_energy)

energy=zeros(nb_structure)
for i=1:nb_structure
    energy[i] = parse(Float64,split(lines[i])[3])
end

file_out=open(string(folder_in,"PIV_1/results",d0,"-",n,"-",m,".dat"),"w")
for i=1:nb_structure
    print("Progress: ",i/nb_structure*100,"%\n")
    for j=i+1:nb_structure
        write(file_out,string(distances[i,j]," ",abs(energy[i]-energy[j])/32*13.605693,"\n"))
    end
end
close(file_out)
