GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

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
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10]
Temperatures=[1750,2000,2500,3000]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10

nb_delta2=5000
nb_space=1000


file_out=open(string(folder_base,"MSD_PV.dat"),"w")
for T in Temperatures
    file_out2=open(string(folder_base,"MSD_PV-",T,"K.dat"),"w")
    for V in Volumes
        file_MSD=string(folder_base,V,"/",T,"K/Data/MSD_coeff.dat")
        file_P=string(folder_base,V,"/",T,"K/Data/Avg_Pressure-BootStrap-nboot_1000.dat")
        if ! isfile(file_MSD) || ! isfile(file_P)
            continue
        end
        file_in=open(file_MSD)
        lines=readlines(file_in)
        close(file_in)

        D=parse(Float64,split(lines[1])[1])

        # Reading file
        file_p=open(file_P)
        lines=readlines(file_p);
        close(file_p)

        P=parse(Float64,split(lines[1])[2])

        write(file_out,string(P," ",T," ",D/6,"\n"))
        write(file_out2,string(P," ",D/6,"\n"))
    end
    write(file_out,string("\n"))
    close(file_out2)
end
close(file_out)
