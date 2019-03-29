GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

Volumes=[9.25]
Temperatures=[2500]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10

nb_delta2=10000
nb_space=500


V=8.82
T=3000

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

barycenter_C=zeros(nb_steps,3)
barycenter_O=zeros(nb_steps,3)
barycenter_all=zeros(nb_steps,3)
for step=1:nb_steps
    for carbon=1:nbC
        for i=1:3
            barycenter_C[step,i] += traj[step].positions[carbon,i]
            barycenter_all[step,i] += traj[step].positions[carbon,i]
        end
    end
    for oxygen=1:nbO
        for i=1:3
            barycenter_O[step,i] += traj[step].positions[nbC+oxygen,i]
            barycenter_all[step,i] += traj[step].positions[nbC+oxygen,i]
        end
    end
    for j=1:3
        barycenter_C[step,j] = barycenter_C[step,j]/nbC
        barycenter_O[step,j] = barycenter_O[step,j]/nbO
        barycenter_all[step,j] = barycenter_all[step,j]/(nbC+nbO)
    end
end

MSD_atoms=zeros(nb_delta2)
MSD_variance=zeros(nb_delta2)
count_2=0
for start_cut=1:nb_space:nb_steps-nb_delta2
    MSD_local=zeros(nb_delta2)
    file_out=open(string(folder_out,"MSD_Atoms_Slide_Total_",count_2+1,".dat"),"w")
    for atom=1:nb_atoms
        positions=traj[start_cut].positions[atom,:]
        for step=1:nb_delta2
            dist=0
            for i=1:3
                dist += ( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )*( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )
            end
            MSD_local[step] += dist
        end
    end
    MSD_local /= (nbC+nbO)
    for i=1:nb_delta2
        write(file_out,string(i*0.005," ",MSD_local[i],"\n"))
        MSD_atoms[i] += MSD_local[i]
        MSD_variance[i] += MSD_local[i]*MSD_local[i]
    end
    global count_2 += 1
    close(file_out)
end
MSD_atoms /= count_2
for i=1:nb_delta2
    MSD_variance[i]= MSD_variance[i]/count_2 - MSD_atoms[i]*MSD_atoms[i]
end

file_out=open(string(folder_out,"MSD_Atoms_Slide_Total_Avg.dat"),"w")
for i=1:nb_delta2
    write(file_out,string(i*0.005," ",MSD_atoms[i]," ",sqrt(MSD_variance[i]),"\n"))
end
close(file_out)
