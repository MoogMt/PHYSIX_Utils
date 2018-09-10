include("CPMD.jl")
include("statistics.jl")

importall CPMD
importall statistics

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

#
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2250,2500,3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures
    for V in volumes

        # defining the path of the files
        folder=string(folder_base,V,"/",T,"K/")

        if isfile( string(folder,"STRESS") )

            # Comoutin pressure
            print( string(folder,"STRESS\n"))
            p=CPMD.readPressure( string(folder,"STRESS") , false , 1)

            # Writting data to disk
            file_pressure=open(string(folder,"presure_time.dat"),"w")
            for step=1:size(p)[1]
                write(file_pressure,string( step ," ", p[step], "\n" ))
            end
            close(file_pressure)

            # Writting results for safe keeping
            write(file,string(statistics.simpleAverage(p)/10," ",T," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
end

close(file)
