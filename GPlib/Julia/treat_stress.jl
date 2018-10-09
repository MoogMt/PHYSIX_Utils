include("CPMD.jl")
include("statistics.jl")

importall CPMD
importall statistics

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

N=96
kb=1.380648*10^-23.

V=
T=2500

run_number=1

fs2ps=0.001
stride=1
timestep=1
unit_sim=stride*fs2ps*timestep
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
