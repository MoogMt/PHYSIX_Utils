GPfolder=string("~/PHYSIX_Utils/GPlib/Julia/")
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

# T,V
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82]
Temperatures=[2000,2500,3000]
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

        for nbrun in runs

            # Input folder
            folder_in=string(folder_local,nbrun,"-run/")

            file_input=string(folder_in,"input")
            if ! isfile(file_input)
                break
            end

            if ! isfile( string(folder_in,"FLAG") ) || ! isfile( string(folder_in,"FLAG2") ) || ! isfile( string(folder_in,"FLAG") )
                print("FLAG spotted in ",V," ",T,"K",)
            end

            if nbrun > 1

                traj_current_curr = filexyz.readFileAtomList( folder_local, nbrun, "-run/TRAJEC_db.xyz")
                traj_current_prev = filexyz.readFileAtomList( folder_local, nbrun-1, "-run/TRAJEC_db.xyz")

                if exp_data.computeRMSD(traj[1],traj2[Int(step)]) < cut_off_rmsd
                    print("Faillure to merge simply at ",nbrun,"\n")
                    handle_out = open( string(folder_local,nbrun,"-run/FLAG4"), "w")
                    write(handle_out,string("CHECK"))
                    close(handle_out)
                    merge = false
                    max_nb_run += 1
                    break
                else
                    times = atom_mod.getNbStep( folder_local, nbrun, "-run/TRAJEC_db.xyz" )
                    nb_total += times
                    push!( total_time, times )
                end
            else
                times = atom_mod.getNbStep( folder_local, nbrun, "-run/TRAJEC_db.xyz" )
                nb_total += times
                push!( total_time, times )
                max_nb_run += 1
            end

        end

        if nb_total < target_length && merge
            print("Total length of simulation is not long enough for merge!")
            print("Length: ",nb_total," step Target: ",target_length," step\n")
        else
            # Merging TRAJEC.xyz
            #---------------------------------------------
            traj_final = Vector{ AtomList }(undef, target_length)
            for i=max_nb_run:-1:1

            end
            filexyz.writeXYZ( folder_local, "TRAJEC_db.xyz" )
            traj_final=0
            # Merging FTRAJ
            #---------------------------------------------
            positions_final = zeros(Real, target_length, 3 )
            velocities_final = zeros(Real, target_length, 3 )
            forces_final = zeros(Real, target_length, 3 )
            for i=max_nb_run:-1:1

            end
            cpmd.writeFtraj( string( folder_local, "FTRAJECTORY_db" ) )
            positions_final=0
            velocities_final=0
            forces_final = 0
            #---------------------------------------------
            # Merging ENERGIES
            temp_final = zeros(Real, target_length )
            epot_final = zeros(Real, target_length )
            etot_final = zeros(Real, target_length )
            msd_final = zeros(Real, target_length )
            comp_time = zeros(Real, target_length )
            for i=max_nb_run:-1:1

            end
            cpmd.writeEnergies( string( folder_local, "ENERGIES_db" ) )
            temp_final = 0
            epot_final = 0
            etot_final = 0
            msd_final = 0
            comp_time = 0
            #---------------------------------------------
            # Merging Stress
            stress_tensor = zeros( Real, target_length,3,3)
            for i=max_nb_run:-1:1

            end
            cpmd.writeStress( string( folder_local, "STRESS_db" ) )
            stress_tensor=0
            #---------------------------------------------
            stress_tensor = zeros( Real, target_length,3,3)
            for i=max_nb_run:-1:1

            end
            cpmd.writePressure( string( folder_local, "Pressure_db" ) )
            pressure=0
            #---------------------------------------------
        end
    end
end
