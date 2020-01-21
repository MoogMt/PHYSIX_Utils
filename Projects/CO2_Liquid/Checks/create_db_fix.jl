GPfolder=string("/home/moogmt/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)

# Aims at fixing stuff
using utils
using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using press_stress
using exp_data

Local_Params=string("/home/moogmt/PHYSIX_Utils/Projects/")
push!(LOAD_PATH, Local_Params )

using computerDB

#==============================================================================#
folder_base = utils.determineFolderPath( computerDB.computers_names, computerDB.computers_pathsCO2 )
if folder_base == false
    print("Computer is not known, add it to the database.")
    exit()
end
#==============================================================================#

# Timesteps
standard_timestep = 40
standard_stride   = 5

# T,V
V = 9.0
T = 3000
run_ = 1

timestep_target=conversion.hatime2fs*standard_timestep*standard_stride

folder_target=string(folder_base,V,"/",T,"K/",run_,"-run/")

file_stress_out = string(folder_target, "STRESS_db")
file_pressure_out = string(folder_target, "Pressure_db.dat")
file_traj_out = string(folder_target,"TRAJEC_db.xyz")
file_ftraj_out = string(folder_target,"FTRAJECTORY_db")
file_energy_out = string(folder_target,"ENERGIES_db")

# Input parameters
file_input=string(folder_target,"input")
stride_stress = cpmd.readIntputStrideStress( file_input )
stride_traj   = cpmd.readIntputStrideTraj( file_input )
timestep_sim  =  cpmd.readInputTimestep( file_input )

# Fixing strides computation
n_stress = round( Int, timestep_target/( timestep_sim*stride_stress ) )
n_traj   = round( Int, timestep_target/( timestep_sim*stride_traj ) )
n_energy = round( Int, timestep_target/( timestep_sim ) )
n_ftraj  = round( Int, timestep_target/( timestep_sim*stride_traj ) )

file_energy_in = string(folder_target,"ENERGIES")
file_trajec_in = string(folder_target,"TRAJEC.xyz")
file_stress_in = string(folder_target,"STRESS")
file_ftrajectory_in = string(folder_target,"FTRAJECTORY")

nb_step_stress = getNbStepStress( file_stress_in )
if nb_step_stress == false
    return false
end
nb_step_stress = utils.nbStepStriding( nb_step_stress, n_stress )
#-------------------------------------------
nb_step_ftraj, nb_atoms_ftraj = getNbStepAtomsFtraj( file_ftrajectory_in )
if nb_step_ftraj == false
    return false
end
nb_step_ftraj  = utils.nbStepStriding( nb_step_ftraj , n_ftraj )
#-------------------------------------------
nb_step_energy = getNbStepEnergies( file_energy_in )
if nb_step_energy == false
    return false
end
nb_step_energy = utils.nbStepStriding( nb_step_energy, n_energy )
#-------------------------------------------
nb_step_traj = filexyz.getNbSteps( file_trajec_in )
if nb_step_traj == false
    return false
end
nb_step_traj = utils.nbStepStriding( nb_step_traj, n_traj )

print("traj_step: ",nb_step_traj,"\n")
print("ftraj_step: ",nb_step_ftraj,"\n")
print("energy_step: ",nb_step_energy,"\n")
print("stress_step: ",nb_step_stress,"\n")
target_length = min( nb_step_traj, nb_step_ftraj, nb_step_energy, nb_step_stress )

nb_ignored=0
#---------------------------------------------
stress_tensor = cpmd.readStress( file_stress_in, n_stress, nb_ignored, target_length )
cpmd.writeStress( file_stress, stress_tensor )
utils.writeData( file_pressure, press_stress.computePressure(stress_tensor) )
stress_tensor=[] # Clearing memory
#---------------------------------------------
filexyz.writeXYZ( file_traj, filexyz.readFileAtomList( file_trajec_in, n_traj, nb_ignored, target_length ) )
#---------------------------------------------
positions,velocities,forces=cpmd.readFtraj( file_ftrajectory_in, n_ftraj, nb_ignored, target_length )
cpmd.writeFtraj( file_ftraj, positions, velocities, forces )
positions=[]
velocities=[]
forces=[]
#---------------------------------------------
temp, epot, etot, msd, comp = cpmd.readEnergies( file_energy_in, n_energy, nb_ignored, target_length )
cpmd.writeEnergies( file_energy, temp, epot, etot, msd, comp )
epot=[]
etot=[]
msd=[]
comp=[]
#---------------------------------------------
