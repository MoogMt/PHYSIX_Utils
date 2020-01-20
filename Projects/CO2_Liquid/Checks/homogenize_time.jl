GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
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

ref_stride_traj=5
ref_timestep=0.4837768653171401


for V in Volumes
    for T in Temperatures
        print("Rolling check on ",V," ",T,"K\n")
        folder_local=string(folder_base,V,"/",T,"K/")

        size_sim=0
        for nbrun in runs

            file_traj=string(folder_local,"TRAJEC.xyz")
            if ! isfile(file_traj)
                print("Break Break Break")
                break
            end

            time_step=cpmd.readInputTimestep( file_input )
            stride_stress=cpmd.readIntputStrideStress( file_input )
            stride_traj=cpmd.readIntputStrideTraj(file_input)

            file_stress=string(folder_in,"STRESS")
            stress,test=cpmd.readStress(file_stress)
            if ! test
                break
            end
            pressure=press_stress.computePressure(stress)
            size_pressure=size(pressure)[1]

            file_energy=string(folder_in,"ENERGIES")
            temperature, e_pot,e_tot,msd,comp_time,test=cpmd.readEnergyFile( file_energy )
            if ! test
                break
            end
            size_energy=size(temperature)[1]

            file_traj=string(folder_in,"TRAJEC.xyz")
            traj,test=filexyz.readFastFile(file_traj)
            if ! test
                break
            end
            size_traj=size(traj)[1]

            filexyz.writeXYZ(file_traj_total,traj_corr)
        end



    end
end
