importall pressure
importall statistics

using PyPlot

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

files=["8.82","9.0","9.05","9.1","9.15","9.2","9.25","9.3","9.35","9.375","9.4","9.5","9.8"]
V=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]
P=Vector{Real}(size(files)[1])
dP=Vector{Real}(size(files)[1])
for i=1:size(files)[1]
    local_file=string(folder,files[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P[i]=statistics.simpleAverage(p)
    dP[i]=sqrt(statistics.simpleMoment(p,2))
end

files_ret=["8.82","9.0r","9.1r","9.2r","9.3r","9.4r","9.5r"]
Vr=[8.82,9.0,9.1,9.2,9.3,9.4,9.5]
Pr=Vector{Real}(size(files_ret)[1])
dPr=Vector{Real}(size(files_ret)[1])
for i=1:size(files_ret)[1]
    local_file=string(folder,files_ret[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    Pr[i]=statistics.simpleAverage(p)
    dPr[i]=sqrt(statistics.simpleMoment(p,2))
end

files_plumed=["9.4+PLUMED3"]
V_Sprint=[9.4]
P_Sprint=Vector{Real}(size(files_plumed)[1])
dP_Sprint=Vector{Real}(size(files_plumed)[1])
for i=1:size(files_plumed)[1]
    local_file=string(folder,files_plumed[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_Sprint[i]=0
    count=0
    for j=1:size(p)[1]
        if p[j] > 0.
            P_Sprint[i]+=p[j]
            count += 1
        end
        P_Sprint[i] = P_Sprint[i]/count
    end
end

plot(V,P,"r.")
plot(Vr,Pr,"b.")
plot(V_Sprint,P_Sprint,"g.")
legend(["3000K - PBE - Up","3000K - PBE - Down","3000K - PBE - SPRINT"])
