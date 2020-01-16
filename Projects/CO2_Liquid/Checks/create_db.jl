GPfolder=string("/home/mathieu/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using press_stress
using exp_data
using utils

function buildingDataBase( folder_target::T1, timestep_target::T2 ) where { T1 <: AbstractString, T2 <: Real }

    # Determining target files paths
    file_input=string(folder_target,"input") # Input - to get timestep + strides
    file_stress=string(folder_target,"STRESS")  # STRESS: contains the stress tensor
    file_traj=string(folder_target,"TRAJEC.xyz") # TRAJEC.xyz: MD Trajectory

    # Getting the stride for the STRESS file
    stride_stress = cpmd.readIntputStrideStress( file_input )
    # Getting the stride for the TRAJEC and FTRAJECTORY files
    stride_traj   = cpmd.readIntputStrideTraj(file_input )
    # Getting the timestep of the simulation
    timestep =  cpmd.readInputTimestep( file_input )

    n_stress =

    file_energy=string(folder_in,"ENERGIES")
    temperature, e_pot,e_tot,msd,comp_time,test=cpmd.readEnergyFile( file_energy )
    if ! test
        return false
    end
    size_data_base=size(temperature)[1]
    temperature=temperature[1:stride_real_data:size_data_base]
    e_pot=e_pot[1:stride_real_data:size_data_base]
    e_tot=e_tot[1:stride_real_data:size_data_base]
    msd=msd[1:stride_real_data:size_data_base]
    comp_time=comp_time[1:stride_real_data:size_data_base]

    stress,test=cpmd.readStress(file_stress)
    if ! test
        return false
    end
    pressure=press_stress.computePressure(stress)
    size_pressure=size(pressure)[1]
    pressure=pressure[1:stride_real_stress:size_pressure]
    size_pressure_actual=size(pressure)[1]

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

#==============================================================================#
computers_names = [
"hp-physix5",
"joliotite",
"OtterCentral"
]
computers_pathsCO2=[
"/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/",
"/media/mathieu/Elements/CO2/",
"/media/moogmt/Element/CO2/"
]
folder_path = utils.determineFolderPath( computers_names, computers_pathsCO2 )
if ! folder_path
    print("Computer is not known, add it to the database.")
    exit()
end
#==============================================================================#

# Timesteps
standard_timestep = 40
standard_stride   = 5

# T,V
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82]
Temperatures=[ 2000, 2500, 3000 ]

# Target Timestep
timestep_target=conversion.hatime2fs*standard_timestep*standard_stride

for V in Volumes
    for T in Temperatures
        n_run = 1
        check = true
        while check
            folder_target=string(folder_base,V,"/",T,"K/")
            if ! isfile( isdir(folder_target) )
                check = false
            end
            if ! buildingDB( folder_target, timestep_target )
                print("Problem at "V," ",T,"K ",n_run," !\n")
            end
        end
    end
end
