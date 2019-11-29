GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))

# Thermodynamical values
Volumes=[8.6,8.82,9.0,9.2,9.3,9.35]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

for cut_off in Cut_Off
    for V in Volumes
        for T in Temperatures
            folder_local=string(folder_base,V,"/",T,"K/Data/")
            if isfile(string(folder_local,"molecules_sizes_cutoff-",cut_off,".dat"))
                print("Progress: V-",V," T-",T,"K","\n")
                file=open(string(folder_local,"molecules_sizes_cutoff-",cut_off,".dat"))
                lines=readlines(file)
                close(file)
                max=0
                for i=1:size(lines)[1]
                    max += parse(Float64,split(lines[i])[3])
                end
                file_out=open(string(folder_local,"size_molecules-",cut_off,"-norm.dat"),"w")
                for i=1:size(lines)[1]
                    write(file_out,string(i," ", parse(Float64,split(lines[i])[3])/max*100,"\n"))
                end
                close(file_out)

                file_out_co2=open(string(folder_local,"co2_ratio.dat"),"w")
                write(file_out_co2,string( parse(Float64,split(lines[3])[3])/max, "\n" ) )
                close(file_out_co2)
            end
        end
    end
end
