GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[10.0,9.8,9.5,9.4,9.375]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]


nbC=32
nbO=nbC*2


print("V=",V," T=",T,"\n")

file_out_map=open(string(folder_base,"map_poisson.dat"),"w")
for T in Temperatures
    for V in Volumes

        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"TRAJEC_wrapped.xyz")
        folder_out=string(folder_in,"Data/")

        if  ! isfile(file)
            continue
        end

        print("Reading Trajectory\n")
        traj=filexyz.readFastFile(file)
        cell=cell_mod.Cell_param(V,V,V)

        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]

        cut_off_bond = 1.75
        max_neigh=5

        data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
        nb_types=size(types)[1]
        states, state_matrices, counts = assignDataToStates( data, nb_types, type_atoms )

        delta=200 # 200 steps = 1ps
        d_delta=50
        occurences_nb=[]
        for step_start=delta:d_delta:nb_steps-delta
            print("Progress: ",step_start/(nb_steps-delta)*100,"%\n")
            for carbon=1:nbC
                occurence=0
                for step=step_start:step_start+delta
                    if state_matrices[1][carbon,step] != 1
                        occurence += 1
                        # Looking up next valid step
                        for next=step+1:step_start+delta
                            if state_matrices[1][carbon,next] == 1
                                step = next-1 # it will get incremeted at end of loop
                            end
                        end
                    end
                end
                push!(occurences_nb,occurence)
            end
        end


        file_out=open(string(folder_out,string("poisson-",delta,".dat")),"w")
        for i=1:size(occurences_nb)[1]
            write(file_out,string(i," ",occurences_nb[i],"\n"))
        end
        close(file_out)

        nb_=size(occurences_nb)[1]

        lambda=sum(occurences_nb)/nb_

        file_in_p=open(string(folder_out,"Avg_Pressure-BootStrap-nboot_1000.dat"))
        lines=readlines(file_in_p)
        close(file_in_p)

        P=parse(Float64,split(lines[1])[2])


        write(file_out_map,string(P," ",T," ",lambda,"\n"))

    end
end
close(file_out_map)
