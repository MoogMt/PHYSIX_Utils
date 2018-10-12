# Loading file
include("contactmatrix.jl")

nbC=32
nbO=64

cut_off=1.75
unit=0.005

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

file_map_3=open(string(folder_out,"3-Coordinances_map.dat"),"w")
file_map_4=open(string(folder_out,"4-Coordinances_map.dat"),"w")

for T in Temperatures
    for V in Volumes

        folder_in=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")

        # Dimers
        #====================================================================================#
        file=string("Avg_Coordinances.dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        file_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/Avg_Pressure-BootStrap-nboot_1000.dat"))
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])

        print(V," ",T,"\n")

        file_read=open(string(folder_in,file));
        lines=readlines(file_read);
        close(file_read);

        C3=parse(Float64, split(lines[1])[2] )
        C4=parse(Float64, split(lines[1])[3] )

        write(file_map_3,string(P," ",T," ",C3,"\n"))
        write(file_map_4,string(P," ",T," ",C4,"\n"))

    end

end

close(file_map_3)
close(file_map_4)
