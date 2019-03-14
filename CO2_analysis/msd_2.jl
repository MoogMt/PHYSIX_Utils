GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

#Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
#Temperatures=[1750,2000,2500,3000]

Volumes=[9.3]
Temperatures=[3000]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10

nb_delta2=5000
nb_space=1000


for V in Volumes
    for T in Temperatures
        print(V," ",T,"\n")
        # Folders
        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"TRAJEC.xyz")
        folder_out=string(folder_in,"Data/")

        if ! isfile(file)
            continue
        end

        # Reading xyz
        print("Computing Data\n")
        traj=filexyz.readFastFile(file)
        cell=cell_mod.Cell_param(V,V,V)

        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]

        MSD_C2=zeros(nb_delta2)
        count_2=0
        for start_cut=1:nb_space:nb_steps-nb_delta2
            MSD_local=zeros(nb_delta2)
            for carbon=1:nbC
                positions=traj[start_cut].positions[carbon,:]
                count=1
                for step=1:nb_delta2
                    dist=0
                    for i=1:3
                        dist += (traj[step+start_cut].positions[carbon,i]-traj[start_cut].positions[carbon,i])*(traj[step+start_cut].positions[carbon,i]-traj[start_cut].positions[carbon,i])
                    end
                    dist=sqrt(dist)
                    MSD_local[step] += dist*dist
                    count+=1
                end
            end
            MSD_local /= nbC
            for i=1:nb_delta2
                MSD_C2[i] += MSD_local[i]
            end
            count_2 += 1
        end
        MSD_C2 /= count_2

        file_out=open(string(folder_out,"MSD_C_slide_total.dat"),"w")
        for i=1:nb_delta2
            write(file_out,string(i*0.005," ",MSD_C2[i],"\n"))
        end
        close(file_out)

        MSD_O2=zeros(nb_delta2)
        count_2=0
        for start_cut=1:nb_space:nb_steps-nb_delta2
            MSD_local=zeros(nb_delta2)
            for oxygen=1:nbO
                positions=traj[start_cut].positions[nbC+oxygen,:]
                for step=1:nb_delta2
                    dist=0
                    for i=1:3
                        dist += (traj[step+start_cut].positions[nbC+oxygen,i]-traj[start_cut].positions[nbC+oxygen,i])*(traj[step+start_cut].positions[nbC+oxygen,i]-traj[start_cut].positions[nbC+oxygen,i])
                    end
                    dist=sqrt(dist)
                    MSD_local[step] += dist*dist
                end
            end
            MSD_local /= nbO
            for i=1:nb_delta2
                MSD_O2[i] += MSD_local[i]
            end
            count_2 += 1
        end
        MSD_O2 /= count_2

        file_out=open(string(folder_out,"MSD_O_slide_total.dat"),"w")
        for i=1:nb_delta2
            write(file_out,string(i*0.005," ",MSD_O2[i],"\n"))
        end
        close(file_out)

        file_out=open(string(folder_out,"MSD_CO_slide_total.dat"),"w")
        for i=1:nb_delta2
            write(file_out,string(i*0.005," ",MSD_C2[i]-MSD_O2[i],"\n"))
        end
        close(file_out)

        MSD_atoms=zeros(nb_delta2)
        count_2=0
        for start_cut=1:nb_space:nb_steps-nb_delta2
            MSD_local=zeros(nb_delta2)
            for atom=1:nb_atoms
                positions=traj[start_cut].positions[atom,:]
                for step=1:nb_delta2
                    dist=0
                    for i=1:3
                        dist += (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i])*(traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i])
                    end
                    dist=sqrt(dist)
                    MSD_local[step] += dist*dist
                end
            end
            MSD_local /= nbO
            for i=1:nb_delta2
                MSD_atoms[i] += MSD_local[i]
            end
            count_2 += 1
        end
        MSD_atoms /= count_2

        file_out=open(string(folder_out,"MSD_Atoms_Slide_Total.dat"),"w")
        for i=1:nb_delta2
            write(file_out,string(i*0.005," ",MSD_atoms[i],"\n"))
        end
        close(file_out)

    end
end
