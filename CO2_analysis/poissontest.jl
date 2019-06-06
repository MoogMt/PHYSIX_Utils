GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.6]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]


nbC=32
nbO=nbC*2

V=9.8
T=3000

print("V=",V," T=",T,"\n")

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

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

unit=0.005

lengths=[]
for carbon=1:nbC
    counting=false
    count_=0
    for step=1:nb_steps
        if ! counting
            # Found event
            if state_matrices[1][carbon,step] != 1
                # We move at the end of the chain
                check=false
                for step_2=1:nb_steps
                    if state_matrices[1][carbon,step_2] == 1
                        step=step_2
                        check=true
                        break
                    end
                end
                # If we did not find it, we go to the next carbon
                if ! check
                    step=nb_steps+1
                end
                # And we start counting to the next event
                counting = true
                count_=1
            end
        else
            # Counting steps to next event
            if state_matrices[1][carbon,step] == 1
                count_ += 1
            else
                if count_ > 1
                    push!(lengths,count_*unit)
                end
                counting=false
                count_=0
            end
        end
    end
end

min_value=unit
max_value=unit
for i=1:size(lengths)[1]
    if max_value < lengths[i]
        global max_value = lengths[i]
    end
end
