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

for V in Volumes
    for T in Temperatures
        print("Merging trajectories in ",V," ",T,"K\n")
        folder_local=string(folder_base,V,"/",T,"K/")

        total_nb_step=0
        total_time=[]

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
            end
        end

        if nb_total < target_length && merge
            print("Total length of simulation is not ")
        end

    end
end
