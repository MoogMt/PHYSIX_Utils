# Loading file
include("contactmatrix.jl")

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

Ry2ev=13.6056980659

folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

energy=open(string(folder_out,"energy_map.dat"),"w")
enthalpy=open(string(folder_out,"enthalpy_map.dat"),"w")

for T in Temperatures
    for V in Volumes

        folder_in=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")

        # Dimers
        #====================================================================================#
        file=string("Avg_EKS-BootStrap-nboot_1000.dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        # Read pressure
        file_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/Avg_Pressure-BootStrap-nboot_1000.dat"))
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])

        # Progress printing
        print(V," ",T,"\n")



        file_read=open(string(folder_in,file));
        lines=readlines(file_read);
        close(file_read);

        EKS=parse(Float64, split(lines[1])[2] )

        write(energy,string(P," ",T," ",EKS*Ry2eV/32,"\n"))
        write(enthalpy,string(P," ",T," ",EKS*Ry2eV/32+P*(10.^9)*V*V*V*(10.0^(-30))*(6.241509*10.^18)/32,"\n"))

    end

end

close(energy)
close(enthalpy)
