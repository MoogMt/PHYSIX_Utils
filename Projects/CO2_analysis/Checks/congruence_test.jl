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
using std_analysis

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# T,V
V=9.05
T=3000
runs=[1,2,3,4]

cut_off_rmsd=0.1

folder_local=string(folder_base,V,"/",T,"K/")

for nbrun in runs

    # Input folder
    folder_in=string(folder_local,nbrun,"-run/")

    file_input=string(folder_in,"input")
    if ! isfile(file_input)
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

    if nbrun > 1
        folder_in_2=string(folder_base,V,"/",T,"K/",nbrun-1,"-run/")
        file_traj_2=string(folder_in_2,"TRAJEC.xyz")
        traj2,test=filexyz.readFastFile(file_traj_2)
        size_traj2=size(traj2)[1]
        connexion_step=-1
        for step=size_traj2:-1:1
            if std_analysis.computeRMSD(traj[1],traj2[step]) < cut_off_rmsd
                connexion_step=step
                break
            end
        end
        if size_traj2 - connexion_step > stride_traj
            print("CHECK FOR CUT, V=",V," T=",T,"K runs:",nbrun-1,"-",nbrun,"\n")
            print("connexion step: ",connexion_step,"\n")
        end
    end

end
