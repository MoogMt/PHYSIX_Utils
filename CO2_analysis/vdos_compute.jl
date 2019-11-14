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

using FFTW

x=range(0,stop=2*pi,step=0.001)
y=sin.(x)
z=dct(y)

file_out=open(string("/home/moogmt/test.dat"),"w")
for i=1:size(z)[1]
    Base.write(file_out,string(i," ",z[i],"\n"))
end
close(file_out)


function vdosFromPosition( file_traj, file_out, max_lag_step )

    traj,cell,test=readFastFile(file_traj)

    if ! sucess
        return zeros(1,1), test
    end



    return vdos, success
end

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=64
nb_atoms=nbC+nbO

T=3000
V=9.8

time_step=0.001
unit_sim=0.5
stride_step=5
dt=time_step*stride_step*unit_sim
dx=0.1 #Angstrom to nm

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC.xyz")
folder_out=string(folder_in,"Data/")

traj,test=readFastFile(file)
velocity=cell_mod.velocityFromPosition(traj,dt,dx)

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
