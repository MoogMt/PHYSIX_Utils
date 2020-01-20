GPfolder=string("/home/moogmt/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)

using utils
using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using press_stress
using exp_data

Local_Params=string("/home/moogmt/PHYSIX_Utils/Projects/")
push!(LOAD_PATH, Local_Params )

using computerDB

#==============================================================================#
folder_path = utils.determineFolderPath( computers_names, computers_pathsCO2 )
if folder_path == false
    print("Computer is not known, add it to the database.")
    exit()
end
#==============================================================================#

# Timesteps
standard_timestep = 40
standard_stride   = 5

# T,V
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82]
Temperatures=[ 2000, 2500, 3000 ]

# Target Timestep
timestep_target=conversion.hatime2fs*standard_timestep*standard_stride

for V in Volumes
    for T in Temperatures
        n_run = 1
        check = true
        while check
            folder_target=string(folder_base,V,"/",T,"K/")
            if ! isfile( isdir(folder_target) )
                check = false
            end
            if ! buildingDB( folder_target, timestep_target )
                print("Problem at "V," ",T,"K ",n_run," !\n")
            end
        end
    end
end
