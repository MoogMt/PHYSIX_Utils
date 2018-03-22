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

files_2K=["8.82","9.0","9.05","9.1","9.2","9.3","9.8"]
V_2K=[8.82,9.0,9.05,9.1,9.2,9.3,9.8]
P_2K=Vector{Real}(size(files_2K)[1])
dP_2K=Vector{Real}(size(files_2K)[1])
for i=1:size(files_2K)[1]
    local_file=string(folder,files_2K[i],"/2000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_2K[i]=statistics.simpleAverage(p)
    dP_2K[i]=sqrt(statistics.simpleMoment(p,2))
end

files_25K=["8.82","9.0","9.05","9.1","9.2","9.3","9.8"]
V_25K=[8.82,9.0,9.05,9.1,9.2,9.3,9.8]
P_25K=Vector{Real}(size(files_25K)[1])
dP_25K=Vector{Real}(size(files_25K)[1])
for i=1:size(files_25K)[1]
    local_file=string(folder,files_25K[i],"/2500K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_25K[i]=statistics.simpleAverage(p)
    dP_25K[i]=sqrt(statistics.simpleMoment(p,2))
end

files_03K=["8.82","9.0","9.2"]
V_03K=[8.82,9.0,9.2]
sizeV=size(files_03K)[1]
P_03K=Vector{Real}(sizeV)
dP_03K=Vector{Real}(sizeV)
for i=1:sizeV
    local_file=string(folder,files_03K[i],"/300K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_03K[i]=statistics.simpleAverage(p)
    dP_03K[i]=sqrt(statistics.simpleMoment(p,2))
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

file_PLUMED=["PLUMED/9.4+PLUMED","PLUMED/9.4+PLUMED3"]
V_PLUMED=[8.82,9.0,9.1,9.2,9.3,9.4,9.5]
P_PLUMED=Vector{Real}(size(file_PLUMED)[1])
dP_PLUMED=Vector{Real}(size(file_PLUMED)[1])
for i=1:size(file_PLUMED)[1]
    local_file=string(folder,file_PLUMED[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_PLUMED[i]=statistics.simpleAverage(p)
    dP_PLUMED[i]=sqrt(statistics.simpleMoment(p,2))
end

figure()
plot(V,P,"r.")
plot(V,P,"r-")
plot(Vr,Pr,"b.")
plot(Vr,Pr,"b-")

figure()
plot(V,P,"r.-")
plot([9.4],P_avg,".")
plot(V_25K,P_25K,"c.-")
plot(V_2K,P_2K,"g.-")
plot(V_03K,P_03K,".-")
legend(["3000K","3000K + SPRINT","2500K","2000K","300K"])
xlabel("V (A**3/atom)")
ylabel("P (kBar)")

plot(V_Sprint,P_Sprint,"g.")
legend(["3000K - PBE - Up","3000K - PBE - Down","3000K - PBE - SPRINT"])
