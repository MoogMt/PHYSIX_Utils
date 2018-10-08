include("contactmatrix.jl")

func="PBE-MT"
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2250,2500,3000]
stride=1
unit=0.005

stop_100 = 20000

nb_box=4

nbC = 32
nbO = 64

for V in Volumes

    for T in Temperatures

        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
        file="TRAJEC.xyz"

        if isfile( string(folder,file) )

            #==============================================================#
            traj = filexyz.read( string(folder,file), stride )
            cell=cell_mod.Cell_param( V, V, V )
            nb_steps=size( traj )[1]
            nb_atoms=size( traj[1].names )[1]
            start_step=Int(1)
            #==============================================================#

            for box=1:nb_box
                #==============================================================#
                bary_origin   = zeros(3)
                bary_step   = zeros(3)
                #==============================================================#

                #==============================================================#
                stop_step=Int(start_step+stop_100/nb_box)
                x0_std_block = traj[start_step].positions[:,:]
                #==============================================================#

                # Barycenter computation
                #==============================================================#
                for atom=1:nb_atoms
                    for i=1:3
                        bary_origin[i] += traj[start_step].positions[atom,i]
                    end
                end
                bary_origin   /= nb_atoms
                #==============================================================#

                #==============================================================#
                file_std=open(string(folder,"MSD-",box,".dat"),"w")
                #==============================================================#

                # Step
                #==============================================================#
                for step=start_step:stop_step
                    #==============================================================#
                    print("V = ",V," T = ",T," Progress: ",step/stop_100*100,"% \n")
                    #==============================================================#

                    # Check to avoid error
                    #==============================================================#
                    if step > nb_steps
                        break
                    end
                    #==============================================================#

                    # Compute local barycenter
                    #==============================================================#
                    bary_step   = zeros(3)
                    for atom=1:nb_atoms
                        for i=1:3
                            bary_step[i] += traj[step].positions[atom,i]
                        end
                    end
                    bary_step  /= nb_atoms
                    #==============================================================#

                    # WRAP
                    #==========================================================#
                    msd_local_std = 0
                    for atom=1:nb_atoms
                        for i=1:3
                            msd_local_std+=((traj[step].positions[atom,i]-traj[1].positions[atom,i])-(bary_step[i]-bary_origin[i]))*((traj[step].positions[atom,i]-traj[1].positions[atom,i])-(bary_step[i]-bary_origin[i]))
                        end
                    end
                    #==========================================================#

                    #==========================================================#
                    write( file_std,   string( step*unit," ", msd_local_std/nb_atoms,"\n"   ) )
                    #==========================================================#

                end
                #==========================================================#
                close( file_std )
                #==========================================================#

                start_step = Int(stop_step+1)
            end


        end

    end

end
