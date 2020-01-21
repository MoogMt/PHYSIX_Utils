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
folder_base = utils.determineFolderPath( computerDB.computers_names, computerDB.computers_pathsCO2 )
if folder_base == false
    print("Computer is not known, add it to the database.")
    exit()
end
#==============================================================================#

# Timesteps
standard_timestep = 40
standard_stride   = 5

# T,V
Volumes=[9.375]
Temperatures=[ 3000 ]


target_timestep=conversion.hatime2fs*standard_timestep*standard_stride

# Target Timestep
for V in Volumes
    for T in Temperatures
        check = true
        run=1
        while check
            folder_target=string(folder_base,V,"/",T,"K/",run,"-run/")
            if ! isdir(folder_target)
                check = false
                continue
            else
                file_stress_out = string(folder_target, "STRESS_db")
                file_pressure_out = string(folder_target, "Pressure_db.dat")
                file_traj_out = string(folder_target,"TRAJEC_db.xyz")
                file_ftraj_out = string(folder_target,"FTRAJECTORY_db")
                file_energy_out = string(folder_target,"ENERGIES_db")
                test=cpmd.buildingDataBase( folder_target, file_stress_out, file_pressure_out, file_traj_out, file_ftraj_out, file_energy_out, target_timestep )
                # buildingDataBase contains all necessary error messages, at least in theory...
                if test
                    print("OK for V=",V," T=",T,"K run=",run," \n")
                else
                    print("Problem for V=",V," T=",T,"K run=",run," \n")
                end
                run += 1
            end
        end
    end
end
