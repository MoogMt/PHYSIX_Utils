# Loading file
include("contactmatrix.jl")

nbC=32
nbO=64

cut_off=1.75
unit=0.005

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

file_CC_occ=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/CC-",cut_off,"_map_occurence.dat"),"w")
file_CC_life=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/CC-",cut_off,"_map_lifetime.dat"),"w")

file_CC_occ_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/CC-",cut_off,"_map_occurence_p.dat"),"w")
file_CC_life_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/CC-",cut_off,"_map_lifetime_p.dat"),"w")

file_dimer_occ=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/dimer-",cut_off,"_map_occurence.dat"),"w")
file_dimer_life=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/dimer-",cut_off,"_map_lifetime.dat"),"w")

file_dimer_occ_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/dimer-",cut_off,"_map_occurence_p.dat"),"w")
file_dimer_life_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/dimer-",cut_off,"_map_lifetime_p.dat"),"w")

file_trimer_occ=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/trimer-",cut_off,"_map_occurence.dat"),"w")
file_trimer_life=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/trimer-",cut_off,"_map_lifetime.dat"),"w")

file_trimer_occ_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/trimer-",cut_off,"_map_occurence_p.dat"),"w")
file_trimer_life_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/trimer-",cut_off,"_map_lifetime_p.dat"),"w")

for T in Temperatures
    for V in Volumes

        folder_in=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")
        folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")

        # Dimers
        #====================================================================================#
        file=string("dimer-",cut_off,".dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        file_p=open(string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/Avg_Pressure-BootStrap-nboot_1000.dat"))
        lines=readlines(file_p);
        close(file_p)
        P=parse(Float64,split(lines[1])[2])

        print("P=",P,"\n")

        print(V," ",T,"\n")

        file_read=open(string(folder_in,file));
        lines=readlines(file_read);
        close(file_read);
        nb_lines=size(lines)[1]

        nb_occurence = 0
        used=zeros(nb_lines)
        life_average=0

        lifes=open(string(folder_out,"lifes-dimer-",cut_off,".dat"),"w")
        for i=1:nb_lines-1
            if used[i] > 0
                continue
            end
            life=0
            nb_occurence += 1
            time_last=parse( Int, split(lines[i])[1] )
            C1=parse( Int, split(lines[i])[2] )
            C2=parse( Int, split(lines[i])[3] )
            O1=parse( Int, split(lines[i])[4] )
            O2=parse( Int, split(lines[i])[5] )
            for j=i+1:nb_lines
                if C1 == parse( Int, split(lines[j])[2] ) && C2 == parse( Int, split(lines[j])[3] ) && O1 == parse( Int, split(lines[j])[4] ) && O2 == parse( Int, split(lines[j])[5] )
                    time_new=parse( Int, split(lines[j])[1] )
                    if time_last - time_new < 100
                        time_last = time_new
                        life += 1
                    else
                        # REBOOT
                        life_average =+ life
                        life=1
                        nb_occurence += 1
                    end
                    used[j] = 1
                end
            end
            life_average =+ life
            write(lifes,string(life*unit," ",life," ",C1," ",C2," ",O1," ",O2,"\n"))
        end
        close(lifes)


        occurence_out=open(string(folder_out,"occurences_dimer-",cut_off,".dat"),"w")
        write(occurence_out,string(nb_occurence,"\n"))
        write(file_dimer_occ,string(T," ",V," ",nb_occurence,"\n"))
        write(file_dimer_occ_p,string(P," ",T," ",nb_occurence,"\n"))
        close(occurence_out)


        life_average_out=open(string(folder_out,"life_average_dimer",cut_off,".dat"),"w")
        write(life_average_out,string(life_average*unit,"\n"))
        write(file_dimer_life,string(T," ",V," ",life_average*unit,"\n"))
        write(file_dimer_life_p,string(P," ",T," ",life_average*unit,"\n"))
        close(life_average_out)
        #====================================================================================#

        #CC
        #====================================================================================#
        file=string("CCbond-",V,"-",T,"-",cut_off,".dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        print(V," ",T,"\n")

        file_read=open(string(folder_in,file));
        lines=readlines(file_read);
        close(file_read);
        nb_lines=size(lines)[1]

        nb_occurence = 0
        used=zeros(nb_lines)
        life_average=0

        lifes=open(string(folder_out,"lifes-CC-",cut_off,".dat"),"w")
        for i=1:nb_lines-1
            if used[i] > 0
                continue
            end
            life=0
            nb_occurence += 1
            time_last=parse( Int, split(lines[i])[1] )
            C1=parse( Int, split( lines[i])[2] )
            C2=parse( Int, split( lines[i])[3] )
            for j=i+1:nb_lines
                if C1 == parse( Int, split(lines[j])[2] ) && C2 == parse( Int, split(lines[j])[3] )
                    time_new=parse( Int, split(lines[j])[1] )
                    if time_last - time_new < 100
                        time_last = time_new
                        life += 1
                    else
                        # REBOOT
                        life_average =+ life
                        life=1
                        nb_occurence += 1
                    end
                    used[j] = 1
                end
            end
            life_average =+ life
            write(lifes,string(life*unit," ",life," ",C1," ",C2,"\n"))
        end
        close(lifes)

        occurence_out=open(string(folder_out,"occurences_CC-",cut_off,".dat"),"w")
        write(occurence_out,string(nb_occurence,"\n"))
        write(file_CC_occ,string(T," ",V," ",nb_occurence,"\n"))
        write(file_CC_occ_p,string(P," ",T," ",nb_occurence,"\n"))
        close(occurence_out)

        life_average_out=open(string(folder_out,"life_average_CC",cut_off,".dat"),"w")
        write(life_average_out,string(life_average*unit,"\n"))
        write(file_CC_life,string(T," ",V," ",life_average*unit,"\n"))
        write(file_CC_life_p,string(P," ",T," ",life_average*unit,"\n"))
        close(life_average_out)
        #====================================================================================#

        # Trimer
        #====================================================================================#
        file=string("trimer-method2-",cut_off,".dat")
        if ! isfile( string(folder_in,file) )
            continue
        end

        print(V," ",T,"\n")

        file_read=open(string(folder_in,file));
        lines=readlines(file_read);
        close(file_read);
        nb_lines=size(lines)[1]

        nb_occurence = 0
        used=zeros(nb_lines)

        life_average=0

        lifes=open(string(folder_out,"lifes-trimer-",cut_off,".dat"),"w")
        for i=1:nb_lines-1
            if used[i] > 0
                continue
            end
            life=0
            nb_occurence += 1
            time_last=parse( Int, split(lines[i])[1] )
            C1=parse( Int, split( lines[i])[2] )
            C2=parse( Int, split( lines[i])[3] )
            C3=parse( Int, split( lines[i])[4] )
            O1=parse( Int, split( lines[i])[5] )
            O2=parse( Int, split( lines[i])[6] )
            O3=parse( Int, split( lines[i])[7] )
            for j=i+1:nb_lines
                if C1 == parse( Int, split(lines[j])[2] ) && C2 == parse( Int, split(lines[j])[3] ) && C3 == parse( Int, split(lines[j])[4] ) && O1 == parse( Int, split(lines[j])[5] ) && O2 == parse( Int, split(lines[j])[6] ) && O3 == parse( Int, split(lines[j])[7] )
                    time_new=parse( Int, split(lines[j])[1] )
                    if time_last - time_new < 100
                        time_last = time_new
                        life += 1
                    else
                        # REBOOT
                        life_average =+ life
                        life=1
                        nb_occurence += 1
                    end
                    used[j] = 1
                end
            end
            life_average =+ life
            write(lifes,string(life*unit," ",life," ",C1," ",C2," ",C3," ",O1," ",O2," ",O3,"\n"))
        end
        close(lifes)

        occurence_out=open(string(folder_out,"occurences_trimer-",cut_off,".dat"),"w")
        write(occurence_out,string(nb_occurence,"\n"))
        write(file_trimer_occ,string(T," ",V," ",nb_occurence,"\n"))
        write(file_trimer_occ_p,string(P," ",T," ",nb_occurence,"\n"))
        close(occurence_out)

        life_average_out=open(string(folder_out,"life_average_trimer",cut_off,".dat"),"w")
        write(life_average_out,string(life_average*unit,"\n"))
        write(file_trimer_life,string(T," ",V," ",life_average*unit,"\n"))
        write(file_trimer_life_p,string(P," ",T," ",life_average*unit,"\n"))
        close(life_average_out)

    end
    write(file_dimer_occ,string("\n"))
    write(file_CC_occ,string("\n"))
    write(file_trimer_occ,string("\n"))
    write(file_dimer_life,string("\n"))
    write(file_CC_life,string("\n"))
    write(file_trimer_life,string("\n"))
    write(file_dimer_occ_p,string("\n"))
    write(file_CC_occ_p,string("\n"))
    write(file_trimer_occ_p,string("\n"))
    write(file_dimer_life_p,string("\n"))
    write(file_CC_life_p,string("\n"))
    write(file_trimer_life_p,string("\n"))
end

close(file_dimer_occ)
close(file_trimer_occ)
close(file_CC_occ)
close(file_dimer_life)
close(file_CC_life)
close(file_trimer_life)

close(file_dimer_occ_p)
close(file_trimer_occ_p)
close(file_CC_occ_p)
close(file_dimer_life_p)
close(file_CC_life_p)
close(file_trimer_life_p)
