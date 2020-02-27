GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion
using fftw
using exp_data

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base = "/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base = "/media/mathieu/Elements/CO2/"

Temperatures = [ 1750, 2000, 2500, 3000 ]
Volumes = [ 10.0, 9.8, 9.5, 9.4, 9.375, 9.35, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6]

rmin=0
dr=0.001

for V in Volumes
    for T in Temperatures

        print("V: ",V," T:",T,"K\n")

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC_fdb_wrapped.xyz")

        if ! isfile( file_in )
            continue
        end

        folder_out=string(folder_in,"Data/Exp/")

        if !isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end

        file_out_gr=string(folder_out,"gr.dat")

        rmax=V/2

        gr = exp_data.computeGr( file_in, file_out_gr, V, rmin, rmax, dr )

        file_out_fq=string(folder_out,"fq.dat")

        rho=96/V^3

        q, fq = exp_data.computeFQ( file_out_fq, gr, rmin, rmax, dr, rho )
    end
end
