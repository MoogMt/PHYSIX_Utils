GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"cpmd.jl"))
include(string(GPfolder,"statistics.jl"))
include(string(GPfolder,"contactmatrix.jl"))

func="PBE-MT"

Volumes=[9.2,9.3,9.35]
Temperatures=[2000,2500]

for V in Volumes
    for T in Temperatures

        folder_input=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
        folder_output=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/Data/")

        # Pressure
        #==============================================================================#
        file="P_all.dat"

        if ! isfile(string(folder_input,file))
            continue
        end

        file=open(string(folder_input,file));
        lines=readlines(file);
        close(file);

        max_step=20000
        start_step=2000

        if size(lines)[1] < max_step + start_step
            max_step=size(lines)[1]-start_step
        end

        start_launch=0

        kbar2GPa=0.1

        avg_P = 0
        #if size(lines)[1]-start_step > 0 && size(lines)[1] > max_step+start_step
        pressure_inst=zeros(max_step)
        for i=1:max_step
            if  i > size(lines)[1]
                break
            end
            pressure_inst[i] = parse( Float64, split(lines[start_step+i-1])[3] ) *kbar2GPa
            avg_P += pressure_inst[i]
        end
        avg_P /= max_step

        # "Normal " average pressure
        file_out=open(string(folder_output,"Avg_Pressure.dat"),"w")
        write(file_out,string(V," ",avg_P,"\n"))
        close(file_out)

        # Histogram
        max=pressure_inst[1]
        min=pressure_inst[1]
        for i=1:max_step
            if pressure_inst[i] > max
                max = pressure_inst[i]
            end
            if pressure_inst[i] < min
                min = pressure_inst[i]
            end
        end

        min=40

        nb_box=50
        d_pressure=(max-min)/nb_box
        hist_p=zeros(nb_box)
        total=0
        for i=1:nb_box
            for pressure in pressure_inst
                if pressure > min + i*d_pressure && pressure < min + (i+1)*d_pressure
                    hist_p[i] += 1
                    total += 1
                end
            end
        end
        hist_p /= total
        file_hist_out=open(string(folder_output,"Hist_Pressure_nbox-",nb_box,".dat"),"w")
        for i=1:nb_box
            write(file_hist_out,string(min+(i+0.5)*d_pressure," ",hist_p[i],"\n"))
        end
        close(file_hist_out)

        # BootStrap
        avg_P = 0
        nb_boot = 1000
        for i=1:nb_boot
            avg_p_local=0
            for i=1:max_step
                avg_p_local += pressure_inst[ Int(trunc(rand()*max_step)+1) ]
            end
            avg_P  += avg_p_local/max_step
        end
        avg_P /= nb_boot

        file_out=open(string(folder_output,"Avg_Pressure-BootStrap-nboot_",nb_boot,".dat"),"w")
        write(file_out,string(V," ",avg_P,"\n"))
        close(file_out)

        # Writting pertinent pressure in file
        file_P_go=open(string(folder_input,"Pressure.dat"),"w")
        for i=1:max_step
            write(file_P_go,string(i," ",pressure_inst[i],"\n"))
        end
        close(file_P_go)

        #end
        #==============================================================================#

        # Energy Class
        #==============================================================================#
        file="EClass_all.dat"
        file=open(string(folder_input,file));
        lines=readlines(file);
        close(file);

        max_step=20000
        start_step=2000

        if size(lines)[1] < max_step + start_step
            max_step=size(lines)[1]-start_step
        end

        start_launch=0
        avg_Eclass = 0
        #if size(lines)[1]-start_step > 0 && size(lines)[1] > max_step+start_step
        eclass_inst=zeros(max_step)
        for i=1:max_step
            if  i > size(lines)[1]
                break
            end
            eclass_inst[i] = parse( Float64, split(lines[start_step+i-1])[3] )
            avg_Eclass += eclass_inst[i]
        end
        avg_P /= max_step

        # "Normal " average pressure
        file_out=open( string( folder_output, "Avg_Eclass.dat" ), "w")
        write( file_out, string( V, " ", avg_Eclass, "\n") )
        close( file_out )

        # Histogram
        max = eclass_inst[1]
        min = eclass_inst[1]
        for i=1:max_step
            if eclass_inst[i] > max
                max = eclass_inst[i]
            end
            if eclass_inst[i] < min
                min = eclass_inst[i]
            end
        end

        nb_box=50
        d_eclass=(max-min)/nb_box
        hist_eclass=zeros(nb_box)
        total=0
        for i=1:nb_box
            for eclass in eclass_inst
                if eclass > min + i*d_eclass && eclass < min + (i+1)*d_eclass
                    hist_eclass[i] += 1
                    total += 1
                end
            end
        end
        hist_eclass /= total
        file_hist_out=open(string(folder_output,"Hist_Eclass_nbox-",nb_box,".dat"),"w")
        for i=1:nb_box
            write(file_hist_out,string(min+(i+0.5)*d_eclass," ",hist_eclass[i],"\n"))
        end
        close(file_hist_out)

        # BootStrap
        avg_Eclass = 0
        nb_boot = 1000
        for i=1:nb_boot
            avg_eclass_local=0
            for i=1:max_step
                avg_eclass_local += eclass_inst[ Int(trunc(rand()*max_step)+1) ]
            end
            avg_Eclass  += avg_eclass_local/max_step
        end
        avg_Eclass /= nb_boot

        file_out=open(string(folder_output,"Avg_Eclass-BootStrap-nboot_",nb_boot,".dat"),"w")
        write(file_out,string(V," ",avg_Eclass,"\n"))
        close(file_out)

        # Writting pertinent pressure in file
        file_eclass_go=open(string(folder_input,"Eclass.dat"),"w")
        for i=1:max_step
            write(file_eclass_go,string(i," ",eclass_inst[i],"\n"))
        end
        close(file_eclass_go)
        #end
        #==============================================================================#



        # Temp
        #==============================================================================#
        file="Temp_all.dat"
        file=open(string(folder_input,file));
        lines=readlines(file);
        close(file);

        max_step=20000
        start_step=2000

        if size(lines)[1] < max_step + start_step
            max_step=size(lines)[1]-start_step
        end

        start_launch=0

        avg_Temp = 0
        #if size(lines)[1]-start_step > 0 && size(lines)[1] > max_step+start_step
        temp_inst=zeros(max_step)
        for i=1:max_step
            if  i > size(lines)[1]
                break
            end
            temp_inst[i] = parse( Float64, split(lines[start_step+i-1])[3] )
            avg_Temp += temp_inst[i]
        end
        avg_Temp /= max_step

        # "Normal " average pressure
        file_out=open( string( folder_output, "Avg_Temp.dat" ), "w")
        write( file_out, string( V, " ", avg_Temp, "\n") )
        close( file_out )

        # Histogram
        max = temp_inst[1]
        min = temp_inst[1]
        for i=1:max_step
            if temp_inst[i] > max
                max = temp_inst[i]
            end
            if temp_inst[i] < min
                min = temp_inst[i]
            end
        end

        nb_box=50
        d_temp=(max-min)/nb_box
        hist_temp=zeros(nb_box)
        total=0
        for i=1:nb_box
            for temp in temp_inst
                if temp > min + i*d_temp && temp < min + (i+1)*d_temp
                    hist_temp[i] += 1
                    total += 1
                end
            end
        end
        hist_temp /= total
        file_hist_out=open(string(folder_output,"Hist_Temp_nbox-",nb_box,".dat"),"w")
        for i=1:nb_box
            write(file_hist_out,string(min+(i+0.5)*d_temp," ",hist_temp[i],"\n"))
        end
        close(file_hist_out)

        # BootStrap
        avg_Temp = 0
        nb_boot = 1000
        for i=1:nb_boot
            avg_temp_local=0
            for i=1:max_step
                avg_temp_local += temp_inst[ Int(trunc(rand()*max_step)+1) ]
            end
            avg_Temp  += avg_temp_local/max_step
        end
        avg_Temp /= nb_boot

        file_out=open(string(folder_output,"Avg_Temp-BootStrap-nboot_",nb_boot,".dat"),"w")
        write(file_out,string(V," ",avg_Temp,"\n"))
        close(file_out)

        # Writting pertinent pressure in file
        file_temp_go=open(string(folder_input,"Temperature.dat"),"w")
        for i=1:max_step
            write(file_temp_go,string(i," ",temp_inst[i],"\n"))
        end
        close(file_temp_go)
        #end
        #==============================================================================#

        # EKS Class
        #==============================================================================#
        file="EKS_all.dat"
        file=open(string(folder_input,file));
        lines=readlines(file);
        close(file);

        max_step=20000
        start_step=2000

        if size(lines)[1] < max_step + start_step
            max_step=size(lines)[1]-start_step
        end


        start_launch=0

        avg_Eks = 0
        #if size(lines)[1]-start_step > 0 && size(lines)[1] > max_step+start_step
        eks_inst=zeros(max_step)
        for i=1:max_step
            if  i > size(lines)[1]
                break
            end
            eks_inst[i] = parse( Float64, split(lines[start_step+i-1])[3] )
            avg_Eks += eks_inst[i]
        end
        avg_Eks /= max_step

        # "Normal " average pressure
        file_out=open( string( folder_output, "Avg_EKS.dat" ), "w")
        write( file_out, string( V, " ", avg_Eks, "\n") )
        close( file_out )

        # Histogram
        max = eks_inst[1]
        min = eks_inst[1]
        for i=1:max_step
            if eks_inst[i] > max
                max = eks_inst[i]
            end
            if eks_inst[i] < min
                min = eks_inst[i]
            end
        end

        nb_box=50
        d_eks=(max-min)/nb_box
        hist_eks=zeros(nb_box)
        total=0
        for i=1:nb_box
            for eks in eks_inst
                if eks > min + i*d_eks && eks < min + (i+1)*d_eks
                    hist_eks[i] += 1
                    total += 1
                end
            end
        end
        hist_eks /= total
        file_hist_out=open(string(folder_output,"Hist_EKS_nbox-",nb_box,".dat"),"w")
        for i=1:nb_box
            write(file_hist_out,string(min+(i+0.5)*d_eks," ",hist_eks[i],"\n"))
        end
        close(file_hist_out)

        # BootStrap
        avg_EKS = 0
        nb_boot = 1000
        for i=1:nb_boot
            avg_eks_local=0
            for i=1:max_step
                avg_eks_local += eks_inst[ Int(trunc(rand()*max_step)+1) ]
            end
            avg_EKS  += avg_eks_local/max_step
        end
        avg_EKS /= nb_boot

        file_out=open(string(folder_output,"Avg_EKS-BootStrap-nboot_",nb_boot,".dat"),"w")
        write(file_out,string(V," ",avg_EKS,"\n"))
        close(file_out)

        # Writting pertinent pressure in file
        file_eks_go=open(string(folder_input,"EKS.dat"),"w")
        for i=1:max_step
            write(file_eks_go,string(i," ",eks_inst[i],"\n"))
        end
        close(file_eks_go)
        #end
        #==============================================================================#

        # Time
        #==============================================================================#
        file="Time_all.dat"
        file=open(string(folder_input,file));
        lines=readlines(file);
        close(file);

        max_step=20000
        start_step=2000

        if size(lines)[1] < max_step + start_step
            max_step=size(lines)[1]-start_step
        end


        start_launch=0

        avg_Time = 0
        time_inst=zeros(max_step)
        for i=1:max_step
            if  i > size(lines)[1]
                break
            end
            time_inst[i] = parse( Float64, split(lines[start_step+i-1])[3] )
            avg_Time += time_inst[i]
        end
        avg_Time /= max_step

        # "Normal " average pressure
        file_out=open( string( folder_output, "Avg_Time.dat" ), "w")
        write( file_out, string( V, " ", avg_Time, "\n") )
        close( file_out )

        # Histogram
        max = time_inst[1]
        min = time_inst[1]
        for i=1:max_step
            if time_inst[i] > max
                max = time_inst[i]
            end
            if time_inst[i] < min
                min = time_inst[i]
            end
        end

        nb_box=50
        d_time=(max-min)/nb_box
        hist_time=zeros(nb_box)
        total=0
        for i=1:nb_box
            for time in time_inst
                if time > min + i*d_time && time < min + (i+1)*d_time
                    hist_time[i] += 1
                    total += 1
                end
            end
        end
        hist_time /= total
        file_hist_out=open(string(folder_output,"Hist_Time_nbox-",nb_box,".dat"),"w")
        for i=1:nb_box
            write(file_hist_out,string(min+(i+0.5)*d_time," ",hist_time[i],"\n"))
        end
        close(file_hist_out)

        # BootStrap
        avg_Time = 0
        nb_boot = 1000
        for i=1:nb_boot
            avg_time_local=0
            for i=1:max_step
                avg_time_local += time_inst[ Int(trunc(rand()*max_step)+1) ]
            end
            avg_Time  += avg_time_local/max_step
        end
        avg_Time /= nb_boot

        file_out=open(string(folder_output,"Avg_Time-BootStrap-nboot_",nb_boot,".dat"),"w")
        write(file_out,string(V," ",avg_Time,"\n"))
        close(file_out)

        # Writting pertinent pressure in file
        file_time_go=open(string(folder_input,"Time.dat"),"w")
        for i=1:max_step
            write(file_time_go,string(i," ",time_inst[i],"\n"))
        end
        close(file_time_go)
        #end
        #==============================================================================#

    end
end
