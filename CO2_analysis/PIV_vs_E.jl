GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"


T=3000
V=8.82

n_run=2

d0=1.75
n=6

folder_in=string(folder_base,V,"/",T,"K/",n_run,"-run/")

print("Computing Data\n")
traj=filexyz.readFastFile(string(folder_in,"TRAJEC.xyz"))
cell=cell_mod.Cell_param(V,V,V)

traj=traj[1:200]
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
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,carbon2),d0,n)
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
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,nbC+oxygen,nbC+oxygen2),d0,n)
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
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,nbC+oxygen),d0,n)
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

file_out=open(string("/home/moogmt/PIV-",d0,"-",n,"-",m,".dat"),"w")
for i=1:size(piv)[1]
    print("Progress : ",i/size(piv)[1]*100,"%\n")
    write(file_out,string(i," "))
    for j=1:size(piv)[2]
        write(file_out,string(piv[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

distances=zeros(nb_structure,nb_structure)
for step1=1:nb_structure
    print("Computing distance matrix: ",step1/nb_structure*100,"%\n")
    for step2=step1+1:nb_structure
        for i=1:size(piv)[1]
            distances[step1,step2] += (piv[i,step1]-piv[i,step2])*(piv[i,step1]-piv[i,step2])
        end
        distances[step1,step2] = sqrt(distances[step1,step2])
        distances[step2,step1] = distances[step1,step2]
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
