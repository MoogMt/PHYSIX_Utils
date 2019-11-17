GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Extract and write information from CPMD energy file

using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using cpmd

func="PBE-MT"

Volumes=[8.82]
Temperatures=[2000]

run_nb=1

fs2ps=0.001
sim_stride=5
unit_analysis=5
unit_sim=sim_stride*fs2ps*unit_analysis

folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/",run_nb,"-run/")
folder_base=string("/home/moogmt/Data/CO2/CO2_AIMD/")

for V in Volumes
    for T in Temperatures

        folder_in=string(folder_base,V,"/",T,"K/",run_nb,"-run/")
        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"ENERGIES")

        if ! isfile(file)
            continue
        end

        temperature, e_ks, e_class, msd, time_data = cpmd.readEnergy( string(folder,file) )

        temp_file=open(string(folder,"Temp_base"),"w")
        eks_file=open(string(folder,"EKS_base"),"w")
        eclass_file=open(string(folder,"EClass_base"),"w")
        msd_file=open(string(folder,"MSD_base"),"w")
        time_file=open(string(folder,"Time_base"),"w")
        for i=1:time_stride:size(temperature)[1]
            write(temp_file,string(i," ",i*unit_target," ",temperature[i],"\n"))
            write(eclass_file,string(i," ",i*unit_target," ",e_class[i],"\n"))
            write(msd_file,string(i," ",i*unit_target," ",msd[i],"\n"))
            write(time_file,string(i," ",i*unit_target," ",time_data[i],"\n"))
            write(eks_file,string(i," ",i*unit_target," ",e_ks[i],"\n"))
        end
        close(temp_file)
        close(eks_file)
        close(eclass_file)
        close(msd_file)
        close(time_file)

    end
end
