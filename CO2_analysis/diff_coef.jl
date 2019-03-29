GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"clustering.jl"))

using LsqFit

@. model(x, p) = p[1]*x

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]
Temperatures=[2000,2500,3000]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10

nb_delta2=5000
nb_space=100

for T in Temperatures
    file_out=open(string(folder_base,"D_",T,".dat"),"w")
    for V in Volumes

        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"TRAJEC.xyz")
        folder_out=string(folder_in,"Data/")

        file_P=string(folder_base,V,"/",T,"K/Data/Avg_Pressure-BootStrap-nboot_1000.dat")

        if ! isfile(file) || ! isfile(file_P)
            continue
        end

        file_p=open(file_P)
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])


        print(V," ",T,"\n")
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

        MSD_local=zeros(nb_delta2)
        count_2=0

        D_avg=0
        D_var=0

        for start_cut=1:nb_space:nb_steps-nb_delta2
            MSD_local=zeros(nb_delta2)
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

            times_MSD=clustering.simpleSequence(nb_delta2)*0.005

            fit = curve_fit(model, times_MSD , MSD_local , [0.4] )
            D_loc=coef(fit)[1]
            D_avg += D_loc
            D_var += D_loc*D_loc
            count_2 += 1
        end

        D_avg=D_avg/count_2
        D_var=D_var/count_2-D_avg*D_avg
        write(file_out,string(P," ",D_avg/6," ",sqrt(D_var),"\n"))

    end
    close(file_out)
end
