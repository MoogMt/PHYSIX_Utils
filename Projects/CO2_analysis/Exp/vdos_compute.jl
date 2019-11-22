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

Volumes=[9.8]
Temperatures=[3000]

cpmd_stride=40
time_step=0.001*conversion.hatime2fs*cpmd_stride
stride_sim=5
dt=time_step*stride_sim

max_lag_frac=0.5
nb_windows=10


for V in Volumes
    for T in Temperatures

        print("V: ",V," T:",T,"K\n")

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC.xyz")
        
        folder_out=string(folder_in,"Data/")
        file_out=string(folder_out,"vdos-frac_",max_lag_frac,"-nbwin_",nb_windows,".dat")

        freq,vdos,test=exp_data.vdosFromPosition( file_in, file_out, max_lag_frac, dt, nb_windows )

        print(size(vdos)[1],"\n")
    end
end
