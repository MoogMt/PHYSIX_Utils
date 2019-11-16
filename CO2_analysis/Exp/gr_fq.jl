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
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

Temperatures=[3000]
Volumes=[8.82]

rmin=0
dr=0.001

for V in Volumes
    for T in Temperatures

        print("V: ",V," T:",T,"K\n")

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC_wrapped.xyz")
        folder_out=string(folder_in,"Data/")

        file_out_gr=string(folder_out,"gr.dat")
        
        rmax=V/2

        gr,test=exp_data.computeGr( file_in, file_out_gr, V, rmin, rmax, dr )

        file_out_fq=string(folder_out,"fq.dat")

        rho=96/V^3

        q,fq=exp_data.computeFQ(file_out_fq,gr,rmin,rmax,dr,rho)

    end
end
