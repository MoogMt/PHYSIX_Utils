GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion
using statistics
using Statistics
using Bootstrap


folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-Godecker/8.82/3000K/")

file_in=open(string(folder_base,"presure_time.dat"))
lines=readlines(file_in)
close(file_in)

nb_point=size(lines)[1]

pressure_gth=zeros(nb_point)
for i=1:nb_point
    pressure_gth[i] = parse(Float64,split(lines[i])[2])/10
end

# Going from 0.5fs timestep to 1fs
pressure_gth=pressure_gth[4000:5:nb_point]
nb_point=size(pressure_gth)[1]
avg_p_gth=0
for i=1:nb_point
    global avg_p_gth += pressure_gth[i]
end
avg_p_gth /= nb_point
print("average pressure: ",avg_p_gth)

# Bootstrap
pressure_avg=0
nb_boot=200
pressure_boot=zeros(nb_boot)
for i=1:nb_boot
    pressure=0
    for j=1:nb_point
        pressure += pressure_gth[ Int(trunc(rand()*nb_point+1)) ] # Conversion kBar to GPa
    end
    pressure_boot[i] = pressure/nb_point
    global pressure_avg += pressure_boot[i]
end
pressure_avg /= nb_boot
pressure_sig = 0
for i=1:nb_boot
    global pressure_sig += (pressure_boot[i]-pressure_avg)^2
end
pressure_sig /= (nb_boot-1)

print("bootstrap pressure: ",pressure_avg,"\n")
print("bootstrap variance: ",pressure_sig,"\n") # weirdly small variance

nb_point=size(pressure_gth)[1]

sizes,avg,var=statistics.blockAverage(pressure_gth,100,2000,100)
file_out=open(string(folder_base,"block_average_pressure.dat"),"w")
for i=1:size(sizes)[1]
    Base.write(file_out,string(sizes[i]," ",avg[i]," ",var[i],"\n"))
end
close(file_out)

n_boot = 2000


#===================================#.

folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/")

file_in=open(string(folder_base,"Pressure.dat"))
lines=readlines(file_in)
close(file_in)

nb_point=size(lines)[1]

pressure_mt=zeros(nb_point)

avg_p_mt=0
for i=1:nb_point
    pressure_mt[i] = parse(Float64,split(lines[i])[2])
    global avg_p_mt += pressure_mt[i]
end
avg_p_mt /= nb_point
print("average pressure: ",avg_p_mt)
n_boot = 2000

bs1 = bootstrap(Statistics.mean, pressure_mt, BalancedSampling(n_boot))
cil = 0.95;
bci1 = confint(bs1, BasicConfInt(cil));
bci2 = confint(bs1, PercentileConfInt(cil));
