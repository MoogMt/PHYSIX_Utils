GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Computes VDOS of a trajectory

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov


function vdosFromPosition( file_traj, file_out, max_lag_step )

    traj,cell,test=readFastFile(file_traj)

    if ! sucess
        return zeros(1,1), test
    end



    return vdos, success
end

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=64
nb_atoms=nbC+nbO

T=3000
V=9.8

time_step=0.001
unit_sim=0.5
dt=time_step*stride_step*unit_sim

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC.xyz")
folder_out=string(folder_in,"Data/")

print("Reading Trajectory\n")


nb_step=size(traj)[1]

nb_step_velo=nb_step-1

# compute the velocities
velocities=zeros(nb_step_velo,nb_atoms,3)
for step=1:nb_step-1
    for atom=1:nb_atoms
        for i=1:3
            velocities[step,atom,i]=(traj[step].positions[atom,i]-traj[step+1].positions[atom,i])/dt
        end
    end
end

# Scalar product des vitesses
velo_scal=zeros(nb_step_velo,nb_atoms)
for atom=1:nb_atoms
    for step=1:nb_step_velo
        for i=1:3
            velo_scal[step,atom] += velocities[step,atom,i]*velocities[step,atom,i]
        end
    end
end

max_lag=5000
nb_step_velo=nb_step-1
# Compute the autocorrelation
autocor_v=zeros(nb_step-1,3)
for atom=1:nb_atoms
    print("Progress: ",atom/nb_atoms*100," %\n")
    autocor_loc=zeros(nb_step_velo)
    # Loop over tau
    for step_lag=1:max_lag
        for step=1:nb_step_velo-step_lag
            autocor_loc[step_lag] += velo_scal[step,atom]*velo_scal[step+step_lag,atom]
        end
        autocor_loc[step_lag] /= (nb_step-step_lag)
    end
    for i=1:nb_step_velo
        autocor_v[i] += autocor_loc[i]
    end
end
autocor_v /= nb_atoms

file_out=open(string(folder_base,"VDOS_test.dat"),"w")
for tau=1:max_lag
    write(file_out,string(tau," ",autocor_v[tau],"\n"))
end
close(file_out)

using GR
using FFTW

autocor_v=autocor_v[1:max_lag]

GR.xlabel("test")

GR.plot(autocor_v[1:500],xlim=(1,500))


GR.plot(dct(autocor_v[1:500])[1:20],"r-")
