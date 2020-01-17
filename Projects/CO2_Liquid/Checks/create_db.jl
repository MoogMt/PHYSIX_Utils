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

function strideData!( data::Vector{T1}, stride::T2 ) where { T1 <: Real, T2 <: Int }
    if stride < 0 || stride > size(data)[0]
        return false
    return data[1:stride:size(data)[0]]

function buildingDataBase( folder_target::T1, timestep_target::T2 ) where { T1 <: AbstractString, T2 <: Real }

    # Determining target files paths
    #---------------------------------------------------------------------------
    file_input=string(folder_target,"input") # Input - to get timestep + strides
    file_stress=string(folder_target,"STRESS")  # STRESS: contains the stress tensor
    file_traj=string(folder_target,"TRAJEC.xyz") # TRAJEC.xyz: MD Trajectory
    #---------------------------------------------------------------------------

    # Reading input file
    #---------------------------------------------------------------------------
    # Getting the stride for the STRESS file
    stride_stress = cpmd.readIntputStrideStress( file_input )
    # Getting the stride for the TRAJEC and FTRAJECTORY files
    stride_traj   = cpmd.readIntputStrideTraj(file_input )
    # Getting the timestep of the simulation
    timestep_sim  =  cpmd.readInputTimestep( file_input )
    #---------------------------------------------------------------------------

    # Computing
    #---------------------------------------------------------------------------
    n_stress = round( Int, timestep_target/( timestep_sim*stride_stress ) )
    n_traj   = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    n_energy = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    n_ftraj  = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    #---------------------------------------------------------------------------

    # Read ENERGIES File
    #---------------------------------------------------------------------------
    file_energy=string(folder_target,"ENERGIES")
    temperature, e_pot, e_tot, msd, comp_time, test = cpmd.readEnergyFile( file_energy )
    if ! test
        return false
    end
    #---------------------------------------------------------------------------

    # Treating ENERGIES data
    #---------------------------------------------------------------------------
    size_data_base=size(temperature)[1]
    temperature=temperature[1:n_energy:size_data_base]
    e_pot=e_pot[1:n_energy:size_data_base]
    e_tot=e_tot[1:n_energy:size_data_base]
    msd=msd[1:n_energy:size_data_base]
    comp_time=comp_time[1:n_energy:size_data_base]
    #---------------------------------------------------------------------------

    # READ STRESS file
    #---------------------------------------------------------------------------
    stress,test=cpmd.readStress(file_stress)
    if ! test
        return false
    end
    #---------------------------------------------------------------------------

    # Treating STRESS data -> Compute P and Stride
    #---------------------------------------------------------------------------
    pressure=press_stress.computePressure(stress)
    size_pressure=size(pressure)[1]
    pressure=pressure[1:n_stress:size_pressure]
    #---------------------------------------------------------------------------

    # Reading TRAJEC.xyz file
    #---------------------------------------------------------------------------
    file_traj=string(folder_in,"TRAJEC.xyz")
    traj,test=filexyz.readFastFile(file_traj)
    if ! test
        return false
    end
    #--------------------------------------------------------------------------

    # Treating TRAJ data
    #--------------------------------------------------------------------------
    size_traj=size(traj)[0]
    traj = traj[1:n_traj:size_traj]
    #--------------------------------------------------------------------------

    return temperature, e_potential, e_total, msd, comp_time, pressure, traj
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
