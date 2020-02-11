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
# Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82]
# Temperatures=[2000,2500,3000]
Volumes = [8.82]
Temperatures = [3000]
runs=[1,2,3,4]

cut_off_rmsd=0.2
target_length=20000

for V in Volumes
    for T in Temperatures
        print("Merging trajectories in ",V," ",T,"K\n")
        folder_local=string(folder_base,V,"/",T,"K/")

        total_nb_step=0
        total_time=[]
        max_nb_run=1
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

        if total_nb_step < target_length && ! merge
            if ! merge
                print("Issues with one of the runs. \n")
            else
                print("Total length of simulation is not long enough for merge! \n")
                print("Length: ",total_nb_step," step Target: ",target_length," step\n")
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
                        traj_final[ target_length-total_time[nb_times]:target_length ] = traj[1:total_time[i]]
                    else
                        traj_final[1:target_length] = traj[total_time[i]+1-target_length:total_time[i] ]
                        check=false
                        break
                    end
                elseif target_length-sum(total_time[i:nb_times]) > 1
                    traj_final[target_length-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = traj[1:total_time[i]]
                else
                    traj_final[1:target_length-sum(total_time[i+1:nb_times])] = traj[total_time[i]-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                    break
                end
            end
            filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb.xyz" ), traj_final )
            cell = cell_mod.Cell_param(V,V,V)
            filexyz.writeXYZ( string( folder_local, "TRAJEC_fdb_wrapped.xyz" ), cell_mod.wrap( traj_final, cell ) )
            # Clean up
            traj_final=0
            # Merging FTRAJ
            #---------------------------------------------
            positions_final = zeros(Real, target_length, 3 )
            velocities_final = zeros(Real, target_length, 3 )
            forces_final = zeros(Real, target_length, 3 )
            check=false
            for i=max_nb_run:-1:1
                positions, velocities, forces = cpmd.readFtraj( string( folder_local, i, "-run/FTRAJECTORY_db" )  )
                if ! check
                    break
                end
                if i == max_nb_run
                    if target_length-total_time[nb_times] > 1
                        positions_final[ target_length-total_time[nb_times]:target_length ] = positions[1:total_time[i]]
                        velocities_final[ target_length-total_time[nb_times]:target_length ] = velocities[1:total_time[i]]
                        forces_final[ target_length-total_time[nb_times]:target_length ] = forces[1:total_time[i]]
                    else
                        positions_final[1:target_length] = positions[total_time[i]+1-target_length:total_time[i] ]
                        velocities_final[1:target_length] = velocities[total_time[i]+1-target_length:total_time[i] ]
                        forces_final[1:target_length] = forces[total_time[i]+1-target_length:total_time[i] ]
                        check=false
                        break
                    end
                elseif target_length-sum(total_time[i:nb_times]) > 1
                    positions_final[target_length-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = positions[1:total_time[i]]
                    velocities_final[target_length-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = velocities[1:total_time[i]]
                    forces_final[target_length-sum(total_time[i:nb_times]):target_length-sum(total_time[i+1:nb_times])] = forces[1:total_time[i]]
                else
                    positions_final[1:target_length-sum(total_time[i+1:nb_times])] = positions[total_time[i]-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                    velocities_final[1:target_length-sum(total_time[i+1:nb_times])] = velocities[total_time[i]-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
                    forces_final[1:target_length-sum(total_time[i+1:nb_times])] = forces[total_time[i]-(target_length-sum(total_time[i+1:nb_times])):total_time[i]]
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
            check=0
            remain=total_nb_step
            for i=max_nb_run:-1:1
                if remain < 1
                    break
                end
                temp,epot,etot,msd,comp = cpmd.readEnergies( string( folder_local, i, "-run/ENERGIES_db" ) )
                if remain > total_time[i]
                    temp_final[target_length-total_time[i]-check:target_length-check] = temp[1:total_time[i]]
                    epot_final[target_length-total_time[i]-check:target_length-check] = epot[1:total_time[i]]
                    etot_final[target_length-total_time[i]-check:target_length-check] = etot[1:total_time[i]]
                    msd_final[target_length-total_time[i]-check:target_length-check] = msd[1:total_time[i]]
                    comp_final[target_length-total_time[i]-check:target_length-check] = comp[1:total_time[i]]
                    remain -= total_time[i]
                else
                    temp_final[target_length-total_time[i]-check:target_length-check] = temp[total_time[i]-remain:total_time[i]]
                    epot_final[target_length-total_time[i]-check:target_length-check] = epot[total_time[i]-remain:total_time[i]]
                    etot_final[target_length-total_time[i]-check:target_length-check] = etot[total_time[i]-remain:total_time[i]]
                    msd_final[target_length-total_time[i]-check:target_length-check] = msd[total_time[i]-remain:total_time[i]]
                    comp_final[target_length-total_time[i]-check:target_length-check] = comp[total_time[i]-remain:total_time[i]]
                    remain = 0
                end
                check += total_time[i]
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
            check=0
            remain=total_nb_step
            for i=max_nb_run:-1:1
                if remain < 1
                    break
                end
                stress_tensor = cpmd.readStress( string( folder_local, i, "-run/STRESS_db" )  )
                if remain > total_time[i]
                    stress_tensor_final[target_length-total_time[i]-check:target_length-check,:,:] = stress_tensor[1:total_time[i],:,:]
                    remain -= total_time[i]
                else
                    stress_tensor_final[target_length-total_time[i]-check:target_length-check,:,:] = stress_tensor[total_time[i]-remain:total_time[i],:,:]
                    remain = 0
                end
                check += total_time[i]
            end
            cpmd.writeStress( string( folder_local, "STRESS_fdb" ), stress_tensor_final )
            #---------------------------------------------
            pressure_final = press_stress.computePressure( stress_tensor_final )
            cpmd.writePressure( string( folder_local, "Pressure_fdb" ), pressure_final )
            # Clean up
            stress_tensor_final=0
            pressure_final=0
            #---------------------------------------------
        end
    end
end
