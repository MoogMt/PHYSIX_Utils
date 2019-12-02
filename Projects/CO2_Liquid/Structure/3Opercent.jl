GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Counts the percent of C with at least three O neighbors
# For all Temperatures and Volumes
# Needs wrapped trajectory files with thermalization
# steps removed

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz

# Target simulations parameters
Volume=[8.82]
Temperature=[2000]
# Analysis parameter
cut_off=[1.75]

# Time units
fs2ps=0.001  # conversion from fs to ps
stride_sim=5 # stride of the units
time_sim=1   # time step of the simulation
unit=time_sim*fs2ps*stride_sim # actual length of analysis of simulation

# Folder stuff
folder_base=string("/home/moogmt/Data/CO2/CO2_AIMD/")

nbC=32
nbO=64
nb_atoms=nbC+nbO

# Loop over all cut-off to test
for cutoff in cut_off
    # Loop over all temperature
    for T in Temperature
        # Loop over all volumes
        for V in Volume

            # Determining file path
            folder=string(folder_base,V,"/",T,"K/")
            file_target=string(folder,"TRAJEC_wrapped.xyz")

            # Checking if file exists
            if ! isfile( file_target )
                break
            end

            # Reading trajectory
            traj = filexyz.readFastFile(file_target)
            cell  = cell_mod.Cell_param( V, V, V )

            # Nb of steps
            nb_steps=size(traj)[1]

            # Open file
            file_out=open(string(folder,"3CO_count_",V,"_",T,"_",cutoff,".dat"),"w")
            # Loop over step
            for step=1:nb_steps
                print("progress: ",step/nb_steps*100,"%\n")
                count = 0
                for carbon=1:nbC
                    # Compute all distances of target carbon
                    distances=zeros(nb_atoms-nbC)
                    for oxygen=nbC+1:nb_atoms
                        distances[oxygen-nbC]=cell_mod.distance(traj[step],cell,carbon,oxygen)
                    end
                    # Sort distances
                    sort!(distances)
                    # If the third neighbor is bonded, all atoms are bonded
                    if distances[3] < cutoff
                        count+=1
                    end
                end
                write(file_out,string(step*unit," ",count/nbC*100,"\n"))
            end
            close(file_out)

        end
    end
end
