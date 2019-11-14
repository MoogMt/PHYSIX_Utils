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
using fftw
using correlation

function vdosFromPosition( file_traj::T1 , max_lag_frac::T2 , to_nm::T3, dt::T4 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real }

    # Reading Trajectory
    traj,test=readFastFile(file_traj)

    if ! test
        return zeros(1,1), zeros(1,1), test
    end

    # Computing velocities
    velocity=cell_mod.velocityFromPosition(traj,dt,dx)

    nb_atoms=size(velocity)[2]
    nb_step=size(velocity)[1]

    # Compute scalar product
    velo_scal=zeros(nb_step,nb_atoms)
    for atom=1:nb_atoms
        for step=1:nb_step
            for i=1:3
                velo_scal[step,atom] += velocity[step,atom,i]*velocity[step,atom,i]
            end
        end
    end

    max_lag=Int(trunc(nb_step*max_lag_frac))

    # Average correlation
    autocorr_avg=zeros(max_lag)
    for atom=1:nb_atoms
        autocorr_avg += correlation.autocorrNorm( velo_scal[:,atom] , max_lag )
    end
    autocorr_avg /= nb_atoms

    # Fourrier Transform
    freq,vdos = fftw.doFourierTransformShift( autocorr_avg, dt )

    return freq, vdos, test
end

function vdosFromPosition( file_traj::T1 , file_out::T2 , max_lag_frac::T3 , to_nm::T4, dt::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real, T5 <: Real }

    # Reading Trajectory
    traj,test=readFastFile(file_traj)

    if ! test
        return test
    end

    # Computing velocities
    velocity=cell_mod.velocityFromPosition(traj,dt,dx)

    nb_atoms=size(velocity)[2]
    nb_step=size(velocity)[1]

    # Compute scalar product
    velo_scal=zeros(nb_step,nb_atoms)
    for atom=1:nb_atoms
        for step=1:nb_step
            for i=1:3
                velo_scal[step,atom] += velocity[step,atom,i]*velocity[step,atom,i]
            end
        end
    end

    max_lag=Int(trunc(nb_step*max_lag_frac))

    # Average correlation
    autocorr_avg=zeros(max_lag)
    for atom=1:nb_atoms
        autocorr_avg += correlation.autocorrNorm( velo_scal[:,atom] , max_lag )
    end
    autocorr_avg /= nb_atoms

    # Fourrier Transform
    freq,vdos = fftw.doFourierTransformShift( autocorr_avg, dt )

    file_o=open(string(file_out),"w")
    for i=1:size(vdos)[1]
        Base.write(file_o,string(freq[i]," ",vdos[i],"\n"))
    end
    close(file_o)

    return test
end

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

T=3000
V=9.8

time_step=0.001
unit_sim=0.5
stride_step=5
dt=time_step*stride_step*unit_sim
dx=0.1 #Angstrom to nm

folder_in=string(folder_base,V,"/",T,"K/")
file_in=string(folder_in,"TRAJEC.xyz")

folder_out=string(folder_in,"Data/")
file_out=string(folder_out,"vdos.dat")

max_lag_frac=0.8
to_nm=0.1

test=vdosFromPosition( file_in, file_out, max_lag_frac, to_nm, dt )
