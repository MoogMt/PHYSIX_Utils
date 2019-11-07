GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Extract and write information from CPMD energy file

using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using cpmd

include(string(GPfolder,"CO2_analysis/markovCO2.jl"))

func="PBE-MT"

# Folder for data
folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/",func,"/")
folder_base=string("/home/moogmt/CO2/CO2_AIMD/")

# Data analysis
max_neigh=5
delta=200 # 200 steps = 1ps
d_delta=50

# Thermo data
Volumes=[8.82]
Temperatures=[3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=64
nb_atoms=nbC+nbO


# Loop over cut-off
for cut_off in Cut_Off

    # Average over the Phase Diagram
    file_out_map=open(string(folder_base,"map_avg",cut_off,".dat"),"w")

    # Loop over the volumes
    for V in Volumes
        # Loop over the temperatures
        for T in Temperatures
            file_out_temp=open(string(folder_base,"avg_",T,"K.dat"),"w")

            print("V=",V," T=",T,"\n")

            # path to input/output file
            folder_in=string(folder_base,V,"/",T,"K/")
            file=string(folder_in,"TRAJEC_wrapped.xyz")
            folder_out=string(folder_in,"Data/")

            # Check existence of file
            if  ! isfile(file)
                continue
            end

            # Read trajectory
            traj=filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            # Nb of steps of the simulation
            nb_steps=size(traj)[1]

            # Building data matrix
            data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
            nb_types=size(types)[1]

            # Building states
            states, state_matrices, counts = assignDataToStates( data, nb_types, type_atoms )

            # Counting occurences
            occurences_nb=[]
            # Loop over steps
            for step_start=delta:d_delta:nb_steps-delta
                print("Progress: ",step_start/(nb_steps-delta)*100,"%\n")
                # Loop over carbons
                for carbon=1:nbC
                    # counting occurence
                    occurence=0
                    for step=step_start:step_start+delta
                        # IF carbon is not a CO2
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
            var=0
            for i=1:nb_
                var += occurences_nb[i]*occurences_nb[i]
            end
            var = var/nb_ - lambda*lambda

            file_in_p=open(string(folder_out,"Avg_Pressure-BootStrap-nboot_1000.dat"))
            lines=readlines(file_in_p)
            close(file_in_p)

            P=parse(Float64,split(lines[1])[2])

            write(file_out_map,string(P," ",T," ",lambda,"\n"))
            write(file_out_temp,string(P," ",lambda," ",sqrt(var),"\n"))

        end
        close(file_out_temp)
    end
    close(file_out_map)

end
