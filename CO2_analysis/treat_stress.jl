GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"cpmd.jl"))
include(string(GPfolder,"statistics.jl"))

folder="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

N=96
kb=1.380648*10^-23.

V=9.3
T=3000

run_number=1

fs2ps=0.001
time_stride=1
timestep=1
unit_sim=time_stride*fs2ps*timestep
unit_target=0.005

stride_press = Int(unit_target/unit_sim)

local_file = string(folder,V,"/",T,"K/",run_number,"-run/STRESS")
p=CPMD.readPressure( local_file , false , stride_press)

start_step=1

treated=string(folder,V,"/",T,"K/",run_number,"-run/STRESS_copy" )
file_p=open( treated, "w" )
for i=start_step:size(p)[1]
    write( file_p, string( i-(start_step-1), " ", (i-(start_step-1))*unit_target, " ", p[i] ,"\n") )
end
close( file_p )
