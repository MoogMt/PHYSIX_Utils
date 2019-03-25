GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"cpmd.jl"))
include(string(GPfolder,"statistics.jl"))
include(string(GPfolder,"contactmatrix.jl"))


Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

for T in Temperatures

    file_pv=open(string(folder_base,"PV-",T,"K.dat"),"w")

    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/Data/")

        file="Avg_Pressure-BootStrap-nboot_1000.dat"
        if ! isfile( string(folder_in,file) )
            continue
        end

        # Reading file
        file_p=open(string(folder_in,file))
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])

        write(file_pv,string(V," ",V*V*V/96," ",P,"\n"))
    end

    close(file_pv)
end
