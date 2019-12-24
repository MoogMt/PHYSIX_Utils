GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using filexyz
using pdb
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

function buildingDB( folder_target::T1 ) where { T1 <: AbstractString }

    file_input=string(folder_target,"input")
    file_stress=string(folder_target,"STRESS")
    file_traj=string(folder_target,"TRAJEC.xyz")

    time_step=cpmd.readInputTimestep( file_input )
    stride_stress=cpmd.readIntputStrideStress( file_input )
    stride_traj=cpmd.readIntputStrideTraj(file_input)

    stress,test=cpmd.readStress(file_stress)
    if ! test
        return false
    end
    pressure=press_stress.computePressure(stress)
    size_pressure=size(pressure)[1]

    file_energy=string(folder_in,"ENERGIES")
    temperature, e_pot,e_tot,msd,comp_time,test=cpmd.readEnergyFile( file_energy )
    if ! test
        return false
    end
    size_energy=size(temperature)[1]

    file_traj=string(folder_in,"TRAJEC.xyz")
    traj,test=filexyz.readFastFile(file_traj)
    if ! test
        return false
    end
    size_traj=size(traj)[1]

    file_out=open(folder_target,"All_data.dat")
    for i=1:nb_step
        for j=1:nb_data
            write(file_out,)
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    return true
end

for V in Volumes
    for T in Temperatures
        for n_run in runs
            folder_target=string(folder_base,V,"/",T,"K/")
            if ! buildingDB( folder_target )
                print("Problem at "V," ",T,"K ",n_run," !\n")
            end
        end
    end
end
