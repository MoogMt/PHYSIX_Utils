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
using conversion
using exp_data

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base = "/media/mathieu/Elements/CO2/"

Temperatures = [ 1750, 2000, 2500, 3000 ]
Volumes = [ 10.0, 9.8, 9.5, 9.4, 9.375, 9.35, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6]

cpmd_stride=40
time_step=0.001*conversion.hatime2fs*cpmd_stride
stride_sim=5
dt=time_step*stride_sim

max_lag_frac=0.5
windows=[1,2,5,10,20,50,100]

for V in Volumes
    for T in Temperatures

        print("V: ",V," T:",T,"K\n")

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC_fdb.xyz")

        if ! isfile( file_in)
            continue
        end

        folder_out=string(folder_in,"Data/Exp/")

        if ! isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end

        file_out = string( folder_out, "vdos_gen.out")

        freq,vdos = exp_data.vdosFromPosition( file_in, file_out, max_lag_frac, dt )

    end
end
