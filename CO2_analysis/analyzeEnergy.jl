GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"cpmd.jl"))
include(string(GPfolder,"statistics.jl"))
include(string(GPfolder,"contactmatrix.jl"))

func="PBE-MT"

V=9.3
T=3000

run_nb=1
timestep=1

fs2ps=0.001
sim_stride=5
unit_sim=sim_stride*fs2ps*timestep
unit_target=0.005

time_stride = Int(unit_target/unit_sim)

folder=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/",run_nb,"-run/")
file="ENERGIES"

temperature, e_ks, e_class, msd, time_data = CPMD.readEnergy( string(folder,file) )

temp_file=open(string(folder,"Temp_base"),"w")
eks_file=open(string(folder,"EKS_base"),"w")
eclass_file=open(string(folder,"EClass_base"),"w")
msd_file=open(string(folder,"MSD_base"),"w")
time_file=open(string(folder,"Time_base"),"w")
j=1
for i=1:time_stride:size(temperature)[1]
    write(temp_file,string(i," ",i*unit_target," ",temperature[i],"\n"))
    write(eclass_file,string(i," ",i*unit_target," ",e_class[i],"\n"))
    write(msd_file,string(i," ",i*unit_target," ",msd[i],"\n"))
    write(time_file,string(i," ",i*unit_target," ",time_data[i],"\n"))
    write(eks_file,string(i," ",i*unit_target," ",e_ks[i],"\n"))
    i+=1
end
close(temp_file)
close(eks_file)
close(eclass_file)
close(msd_file)
close(time_file)
