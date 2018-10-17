# Loading file
include("contactmatrix.jl")

nbC=32
nbO=64

cut_off=1.75
unit=0.005

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

file_map_2C=open(string(folder_out,"2-CoordinancesC_map.dat"),"w")
file_map_3C=open(string(folder_out,"3-CoordinancesC_map.dat"),"w")
file_map_4C=open(string(folder_out,"4-CoordinancesC_map.dat"),"w")

file_map_1O=open(string(folder_out,"1-CoordinancesO_map.dat"),"w")
file_map_2O=open(string(folder_out,"2-CoordinancesO_map.dat"),"w")

for T in Temperatures
    for V in Volumes

        folder_in=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")

        # Dimers
        #====================================================================================#
        file=string("Avg_CoordinancesC.dat")
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

        write(file_map_2C,string(P," ",T," ",1-C3-C4,"\n"))
        write(file_map_3C,string(P," ",T," ",C3,"\n"))
        write(file_map_4C,string(P," ",T," ",C4,"\n"))

        file=string("Avg_CoordinancesO.dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        file_read2=open(string(folder_in,file));
        lines=readlines(file_read2);
        close(file_read2);

        O1=parse(Float64, split(lines[1])[2] )

        write(file_map_1O,string(P," ",T," ",O1,"\n"))
        write(file_map_2O,string(P," ",T," ",1-O2,"\n"))

    end

end

close(file_map_2C)
close(file_map_3C)
close(file_map_4C)

close(file_map_1O)
close(file_map_2O)
