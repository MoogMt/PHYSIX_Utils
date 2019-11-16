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
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

Temperatures=[2500]
Volumes=[8.82]

time_step=0.001
stride_sim=5
dt=time_step*stride_sim
dx=0.1 #Angstrom to nm

for V in Volumes
    for T in Temperatures

        print("V: ",V," T:",T,"K\n")

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC.xyz")

        max_lag_frac=0.5
        to_nm=1

        folder_out=string(folder_in,"Data/")
        file_out=string(folder_out,"vdos-",max_lag_frac,".dat")

        freq,vdos,test=vdosFromPosition( file_in, file_out, max_lag_frac, to_nm, dt )

    end
end
