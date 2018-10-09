include("contactmatrix.jl")
include("CPMD.jl")
include("statistics.jl")

importall CPMD

func="PBE-MT"

V=10.0
T=3000

run_nb=2
timestep=1

fs2ps=0.001
stride=1
unit_sim=stride*fs2ps*timestep
unit_target=0.005

stride = Int(unit_target/unit_sim)

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/",run_nb,"-run/")
file="ENERGIES"

temperature, e_ks, e_class, msd, time = CPMD.readEnergy( string(folder,file) )

temp_file=open(string(folder,"Temp_base"),"w")
eks_file=open(string(folder,"EKS_base"),"w")
eclass_file=open(string(folder,"EClass_base"),"w")
msd_file=open(string(folder,"MSD_base"),"w")
time_file=open(string(folder,"Time_base"),"w")
j=1
for i=1:stride:size(temperature)[1]
    write(temp_file,string(j," ",j*unit_target," ",temperature[i],"\n"))
    write(eclass_file,string(j," ",j*unit_target," ",e_class[i],"\n"))
    write(msd_file,string(j," ",j*unit_target," ",msd[i],"\n"))
    write(time_file,string(j," ",j*unit_target," ",time[i],"\n"))
    write(eks_file,string(j," ",j*unit_target," ",e_ks[i],"\n"))
    j+=1
end
close(temp_file)
close(eks_file)
close(eclass_file)
close(msd_file)
close(time_file)
