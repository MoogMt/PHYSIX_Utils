include("CPMD.jl")
include("statistics.jl")

importall CPMD
importall statistics

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#
#
volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10]
temperatures=[2000,2250,2500,2750,3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/96," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)

include("CPMD.jl")
include("statistics.jl")

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/Return/"
volumes=[9.0,9.2,9.3,9.4,9.5]
temperatures=[3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/96," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)

include("CPMD.jl")
include("statistics.jl")

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/PLUMED/"
volumes=[8.82,9.05,9.4]
temperatures=[2000,2500,3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/96," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)

include("CPMD.jl")
include("statistics.jl")

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/Super/"
volumes=[13.23,13.9875]
temperatures=[3000]

nbC=108
nbO=2*nbC
nb_atoms=nbC+nbO

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/nb_atoms," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)

include("CPMD.jl")
include("statistics.jl")

importall CPMD
importall statistics

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-Godecker/"
volumes=[8.82,9.0,9.1,9.3,9.35,9.4,9.8]
temperatures=[3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/96," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)


include("CPMD.jl")
include("statistics.jl")

importall CPMD
importall statistics

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/BLYP/"
volumes=[8.82,9.0,9.05,9.1,9.2,9.3,9.35,9.4,9.5,9.8]
temperatures=[2000,2500,3000]

file=open(string(folder_base,"pressure_all.dat"),"w")

for T in temperatures

    file2=open(string(folder_base,"pressure-",T,"K.dat"),"w")
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
            write(file2,string(V*V*V/96," ",statistics.simpleAverage(p)/10," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
            print(string(statistics.simpleAverage(p)," ",sqrt(statistics.simpleMoment(p,2)),"\n"))
        end

    end
    close(file2)
end

close(file)
