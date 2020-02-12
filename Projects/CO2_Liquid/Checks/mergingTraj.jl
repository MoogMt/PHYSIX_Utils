GPfolder=string("/home/mathieu/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using cpmd
using press_stress
using exp_data

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/media/moogmt/Elements/CO2/"
folder_base="/media/mathieu/Elements/CO2/"

# T,V
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperatures=[1750]
runs=[1,2,3,4]

cut_off_rmsd=0.2
target_length=20000

for V in Volumes
    for T in Temperatures
        print("Merging trajectories in ",V," ",T,"K\n")
        folder_local=string(folder_base,V,"/",T,"K/")

        if ! isdir(folder_local)
            continue
        end

        total_nb_step=0
        total_time=[]
        max_nb_run=0
        merge=true

        for nbrun in runs
            # Input folder
            folder_in=string(folder_local,nbrun,"-run/")
            file_input=string(folder_in,"input")
            if ! isfile(file_input)
                break
            end
            if isfile( string(folder_in,"FLAG") ) || isfile( string(folder_in,"FLAG2") ) || isfile( string(folder_in,"FLAG3") )
                print("FLAG spotted in ",V," ",T,"K\n",)
                if isfile( string(folder_in,"FLAG") )
                    print("FLAG\n")
                end
                if isfile( string(folder_in,"FLAG2") )
                    print("FLAG2\n")
                end
                if isfile( string(folder_in,"FLAG3") )
                    print("FLAG3\n")
                end
                merge = false
            end
            if nbrun > 1

                traj_current_curr = filexyz.readFileAtomList( string( folder_local, nbrun, "-run/TRAJEC_db.xyz") )
                traj_current_prev = filexyz.readFileAtomList( string( folder_local, nbrun-1, "-run/TRAJEC_db.xyz") )

                if traj_current_prev == false || traj_current_curr == false
                    break
                end
                rmsd = exp_data.computeRMSD(traj_current_curr[1],traj_current_prev[size(traj_current_prev)[1]] )
                if rmsd > cut_off_rmsd
                    print("Faillure to merge simply at ",nbrun,"\n")
                    print("RMSD ",rmsd,"\n")
                    handle_out = open( string(folder_local,nbrun,"-run/FLAG4"), "w")
                    write(handle_out,string("CHECK"))
                    close(handle_out)
                    merge = false
                    max_nb_run += 1
                    break
                else
                    times = filexyz.getNbSteps( string(folder_local, nbrun, "-run/TRAJEC_db.xyz" ) )
                    total_nb_step += times
                    max_nb_run += 1
                    push!( total_time, times )
                end
            else
                times = filexyz.getNbSteps( string(folder_local, nbrun, "-run/TRAJEC_db.xyz" ) )
                total_nb_step += times
                push!( total_time, times )
                max_nb_run += 1
            end
        end

        nb_times=size(total_time)[1]

        if ! merge
            print("Major issue somewhere here (",V," ",T,"K), moving on.\n")
            continue
        end

        # If only one simulation
        if max_nb_run == 1
            # Shorting TRAJEC.xyz
            #---------------------------------------------
            check=true
            if total_time[1]+1-target_length < 1
                print("Only one trajectory and it is to short to get to the target length.\n")
                print("Target length: ",target_length,"\n")
                print("Size trajectory: ",total_time[1]+1,"\n")
                check=false
                continue
            else
                traj = filexyz.readFileAtomList( string( folder_local, 1, "-run/TRAJEC_db.xyz" ) )
                traj_final = traj[total_time[1]+1-target_length:total_time[1] ]
                filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb.xyz" ), traj_final )
                cell = cell_mod.Cell_param(V,V,V)
                filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb_wrapped.xyz" ), cell_mod.wrap( traj_final,cell ) )
                nb_atoms = size(traj[1].names)[1]
                # Clean up
                traj=0
                # Shorting FTRAJECTORY
                #---------------------------------------------
                positions, velocities, forces = cpmd.readFtraj( string( folder_local, 1, "-run/FTRAJECTORY_db" )  )
                positions_final = positions[total_time[1]+1-target_length:total_time[1],:,:]
                velocities_final = velocities[total_time[1]+1-target_length:total_time[1],:,:]
                forces_final = forces[total_time[1]+1-target_length:total_time[1],:,:]
                cpmd.writeFtraj( string( folder_local, "FTRAJECTORY_fdb" ), positions_final, velocities_final, forces_final )
                # Clean up
                positions=0
                velocities=0
                forces=0
                positions_final=0
                velocities_final=0
                forces_final = 0
                # Shorting ENERGIES
                #---------------------------------------------
                temp,epot,etot,msd,comp = cpmd.readEnergies( string( folder_local, 1, "-run/ENERGIES_db" )  )
                temp_final = temp[total_time[1]+1-target_length:total_time[1] ]
                epot_final = epot[total_time[1]+1-target_length:total_time[1] ]
                etot_final = etot[total_time[1]+1-target_length:total_time[1] ]
                msd_final = msd[total_time[1]+1-target_length:total_time[1] ]
                comp_final = comp[total_time[1]+1-target_length:total_time[1] ]
                cpmd.writeEnergies( string( folder_local, "ENERGIES_fdb" ), temp_final, epot_final, etot_final, msd_final, comp_final )
                temp = 0
                epot = 0
                etot = 0
                msd = 0
                comp = 0
                temp_final = 0
                epot_final = 0
                etot_final = 0
                msd_final = 0
                comp_time = 0
                # Shorting STRESS
                #---------------------------------------------
                stress_tensor = cpmd.readStress( string( folder_local, 1, "-run/STRESS_db" )   )
                stress_tensor_final = stress_tensor[total_time[1]+1-target_length:total_time[1],:,:]
                stress_tensor=0
                cpmd.writeStress( string( folder_local, "STRESS_fdb" ), stress_tensor_final )
                pressure_final = press_stress.computePressure( stress_tensor_final )
                stress_tensor_final = 0
                press_stress.writePressure( string( folder_local, "Pressure_fdb" ), pressure_final )
                pressure_final = 0
            end
        else
            if total_nb_step < target_length && ! merge
                if ! merge
                    print("Issues with one of the runs. \n")
                    break
                else
                    print("Total length of simulation is not long enough for merge! \n")
                    print("Length: ",total_nb_step," step Target: ",target_length," step\n")
                    break
                end
            else
                # Merging TRAJEC.xyz
                #---------------------------------------------
                check=true
                traj_final = Vector{ AtomList }(undef, target_length)
                for i=max_nb_run:-1:1
                    traj = filexyz.readFileAtomList( string( folder_local, i, "-run/TRAJEC_db.xyz" ) )
                    if ! check
                        break
                    end
                    if i == max_nb_run
                        if target_length-total_time[nb_times] > 1
                            traj_final[ target_length-total_time[nb_times]+1:target_length ] = traj[1:total_time[i]]
                        else
                            traj_final[1:target_length] = traj[total_time[i]+1-target_length:total_time[i] ]
                            check=false
                            break
                        end
                    elseif target_length-sum(total_time[i:nb_times]) > 1
                        traj_final[target_length-sum(total_time[i:nb_times])+1:target_length-sum(total_time[i+1:nb_times])] = traj[1:total_time[i]]
                    else
                        traj_final[1:target_length-sum(total_time[i+1:nb_times])] = traj[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        break
                    end
                end
                filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb.xyz" ), traj_final )
                cell = cell_mod.Cell_param(V,V,V)
                filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb_wrapped.xyz" ), cell_mod.wrap( traj_final, cell ) )
                nb_atoms = size(traj_final[1].names)[1]
                # Clean up
                traj_final=0
                # Merging FTRAJ
                #---------------------------------------------
                positions_final = zeros(Real, target_length, nb_atoms, 3 )
                velocities_final = zeros(Real, target_length, nb_atoms, 3 )
                forces_final = zeros(Real, target_length, nb_atoms, 3 )
                check=true
                for i=max_nb_run:-1:1
                    positions, velocities, forces = cpmd.readFtraj( string( folder_local, i, "-run/FTRAJECTORY_db" )  )
                    if ! check
                        break
                    end
                    if i == max_nb_run
                        if target_length-total_time[nb_times] > 1
                            positions_final[ target_length-total_time[nb_times]+1:target_length,:,:] = positions[1:total_time[i],:,:]
                            velocities_final[ target_length-total_time[nb_times]+1:target_length,:,:] = velocities[1:total_time[i],:,:]
                            forces_final[ target_length-total_time[nb_times]+1:target_length,:,:] = forces[1:total_time[i],:,:]
                        else
                            positions_final[1:target_length,:,:] = positions[total_time[i]+1-target_length:total_time[i],:,:]
                            velocities_final[1:target_length,:,:] = velocities[total_time[i]+1-target_length:total_time[i],:,:]
                            forces_final[1:target_length,:,:] = forces[total_time[i]+1-target_length:total_time[i],:,:]
                            check=false
                            break
                        end
                    elseif target_length-sum(total_time[i:nb_times]) > 1
                        positions_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times]),:,:] = positions[1:total_time[i],:,:]
                        velocities_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times]),:,:] = velocities[1:total_time[i],:,:]
                        forces_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times]),:,:] = forces[1:total_time[i],:,:]
                    else
                        positions_final[1:target_length-sum(total_time[i+1:nb_times]),:,:] = positions[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i],:,:]
                        velocities_final[1:target_length-sum(total_time[i+1:nb_times]),:,:] = velocities[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i],:,:]
                        forces_final[1:target_length-sum(total_time[i+1:nb_times]),:,:] = forces[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i],:,:]
                        break
                    end
                end
                cpmd.writeFtraj( string( folder_local, "FTRAJECTORY_fdb" ), positions_final, velocities_final, forces_final )
                # Clean up
                positions_final=0
                velocities_final=0
                forces_final = 0
                #---------------------------------------------
                # Merging ENERGIES
                temp_final = zeros(Real, target_length )
                epot_final = zeros(Real, target_length )
                etot_final = zeros(Real, target_length )
                msd_final = zeros(Real, target_length )
                comp_final = zeros(Real, target_length )
                check=true
                for i=max_nb_run:-1:1
                    temp,epot,etot,msd,comp = cpmd.readEnergies( string( folder_local, i, "-run/ENERGIES_db" )  )
                    if ! check
                        break
                    end
                    if i == max_nb_run
                        if target_length-total_time[nb_times] > 1
                            temp_final[ target_length+1-total_time[nb_times]:target_length ] = temp[1:total_time[i]]
                            epot_final[ target_length+1-total_time[nb_times]:target_length ] = epot[1:total_time[i]]
                            etot_final[ target_length+1-total_time[nb_times]:target_length ] = etot[1:total_time[i]]
                            msd_final[ target_length+1-total_time[nb_times]:target_length ] = msd[1:total_time[i]]
                            comp_final[ target_length+1-total_time[nb_times]:target_length ] = comp[1:total_time[i]]
                        else
                            temp_final[1:target_length] = temp[total_time[i]+1-target_length:total_time[i] ]
                            epot_final[1:target_length] = epot[total_time[i]+1-target_length:total_time[i] ]
                            etot_final[1:target_length] = etot[total_time[i]+1-target_length:total_time[i] ]
                            msd_final[1:target_length] = msd[total_time[i]+1-target_length:total_time[i] ]
                            comp_final[1:target_length] = comp[total_time[i]+1-target_length:total_time[i] ]
                            check=false
                            break
                        end
                    elseif target_length-sum(total_time[i:nb_times]) > 1
                        temp_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = temp[1:total_time[i]]
                        epot_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = epot[1:total_time[i]]
                        etot_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = etot[1:total_time[i]]
                        msd_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = msd[1:total_time[i]]
                        comp_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = comp[1:total_time[i]]
                    else
                        temp_final[1:target_length-sum(total_time[i+1:nb_times])] = temp[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        epot_final[1:target_length-sum(total_time[i+1:nb_times])] = epot[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        etot_final[1:target_length-sum(total_time[i+1:nb_times])] = etot[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        msd_final[1:target_length-sum(total_time[i+1:nb_times])] = msd[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        comp_final[1:target_length-sum(total_time[i+1:nb_times])] = comp[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                        break
                    end
                end
                cpmd.writeEnergies( string( folder_local, "ENERGIES_fdb" ), temp_final, epot_final, etot_final, msd_final, comp_final )
                # Clean up
                temp_final = 0
                epot_final = 0
                etot_final = 0
                msd_final = 0
                comp_time = 0
                #---------------------------------------------
                # Merging Stress / Pressure
                stress_tensor_final = zeros( Real, target_length,3,3)
                check=true
                for i=max_nb_run:-1:1
                    stress_tensor = cpmd.readStress( string( folder_local, i, "-run/STRESS_db" )   )
                    if ! check
                        break
                    end
                    if i == max_nb_run
                        if target_length-total_time[nb_times] > 1
                            stress_tensor_final[target_length+1-total_time[nb_times]:target_length,:,:] = stress_tensor[1:total_time[i],:,:]
                        else
                            stress_tensor_final[1:target_length,:,:] = stress_tensor[total_time[i]+1-target_length:total_time[i],:,:]
                            check=false
                            break
                        end
                    elseif target_length-sum(total_time[i:nb_times]) > 1
                        stress_tensor_final[target_length+1-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times]),:,:] = stress_tensor[1:total_time[i],:,:]
                    else
                        stress_tensor_final[1:target_length-sum(total_time[i+1:nb_times]),:,:] = stress_tensor[total_time[i]+1-(target_length-sum(total_time[i+1:nb_times])):total_time[i],:,:]
                        break
                    end
                end
                cpmd.writeStress( string( folder_local, "STRESS_fdb" ), stress_tensor_final )
                #---------------------------------------------
                pressure_final = press_stress.computePressure( stress_tensor_final )
                press_stress.writePressure( string( folder_local, "Pressure_fdb" ), pressure_final )
                # Clean up
                stress_tensor_final=0
                pressure_final=0
                #---------------------------------------------
            end
        end
    end
end
