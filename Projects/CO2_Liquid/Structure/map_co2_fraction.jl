# Loading file
include("contactmatrix.jl")

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

file=open(string(folder_base,"map_co2_fraction.dat"),"w")

for T in Temperatures
    count=0
    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/Data/")

        file_name="co2_ratio.dat"
        if ! isfile( string(folder_in,file_name) )
            continue
        end

        file_p=open(string(folder_in,"Avg_Pressure-BootStrap-nboot_1000.dat"))
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])

        file_co2=open(string(folder_in,file_name))
        lines=readlines(file_co2)
        close(file_co2)
        co2_ratio=parse(Float64,split(lines[1])[1])

        write(file,string(P," ",T," ",co2_ratio,"\n"))
        count += 1
    end
    if count != 0
        write(file,"\n")
    end
end

close(file)
