# Aims at fixing stuff
using utils
using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using press_stress
using exp_data

folder_base = string("/media/mathieu/Elements/CO2/")

# T,V
Volumes=[ 9.05 ]
Temperatures=[ 3000 ]

# Target Timestep
for V in Volumes
    for T in Temperatures
        run_=1
        check = true
        while check
            folder_target = string(folder_base,V,"/",T,"K/",run_,"-run/")
            if isdir(folder_target)
                if isfile( string(folder_target,"FLAG") )
                    cpmd.relaunchRunTrajec( folder_target, string(folder_target,"input_restart") )
                else
                    print("V: ",V," T: ",T,"K\n")
                    cpmd.relaunchRunFtraj( folder_target, string(folder_target,"input_restart") )
                end
                print("Ok for V:",V," T:",T,"K run:",run_,"\n")
                run_ += 1
            else
                check = false
            end
        end
    end
end
