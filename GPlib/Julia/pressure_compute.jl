importall CPMD
importall statistics

using PyPlot

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

N=96
kb=1.380648*10^-23.
T=3000

files=["8.82","9.0","9.05","9.1","9.15","9.2","9.25","9.3","9.35","9.375","9.4","9.5","9.8"]
V=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]
# for i=1:size(V)[1]
#     V[i]=V[i]*V[i]*V[i]/96
# end
P=Vector{Real}(size(files)[1])
dP=Vector{Real}(size(files)[1])
for i=1:size(files)[1]
    local_file=string(folder,files[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1) c
    P[i]=statistics.simpleAverage(p)+N*kb*T/(V*10^-10.)^3.
    dP[i]=sqrt(statistics.simpleMoment(p,2))
end

files_2K=["8.82","9.0","9.05","9.1","9.2","9.3","9.8"]
V_2K=[8.82,9.0,9.05,9.1,9.2,9.3,9.8]
for i=1:size(V_2K)[1]
    V_2K[i]=V_2K[i]*V_2K[i]*V_2K[i]/96
end
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
for i=1:size(V_25K)[1]
    V_25K[i]=V_25K[i]*V_25K[i]*V_25K[i]/96
end
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
for i=1:size(V_03K)[1]
    V_03K[i]=V_03K[i]*V_03K[i]*V_03K[i]/96
end
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
for i=1:size(Vr)[1]
    Vr[i]=Vr[i]*Vr[i]*Vr[i]/96
end
Pr=Vector{Real}(size(files_ret)[1])
dPr=Vector{Real}(size(files_ret)[1])
for i=1:size(files_ret)[1]
    local_file=string(folder,files_ret[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    Pr[i]=statistics.simpleAverage(p)
    dPr[i]=sqrt(statistics.simpleMoment(p,2))
end

file_PLUMED=["PLUMED/9.4+PLUMED","PLUMED/9.4+PLUMED3"]
V_PLUMED=[9.4,9.4]
for i=1:size(V_PLUMED)[1]
    V_PLUMED[i]=V_PLUMED[i]*V_PLUMED[i]*V_PLUMED[i]/96
end
P_PLUMED=Vector{Real}(size(file_PLUMED)[1])
dP_PLUMED=Vector{Real}(size(file_PLUMED)[1])
for i=1:size(file_PLUMED)[1]
    local_file=string(folder,file_PLUMED[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_PLUMED[i]=statistics.simpleAverage(p)
    dP_PLUMED[i]=sqrt(statistics.simpleMoment(p,2))
end

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/BLYP/"
files=["8.82","9.05","9.2","9.3","9.35","9.4","9.5","9.8"]
V_BLYP=[8.82,9.05,9.2,9.3,9.35,9.4,9.5,9.8]
for i=1:size(V_BLYP)[1]
    V_BLYP[i]=V_BLYP[i]*V_BLYP[i]*V_BLYP[i]/96
end
P_BLYP=Vector{Real}(size(files)[1])
dP_BLYP=Vector{Real}(size(files)[1])
for i=1:size(files)[1]
    local_file=string(folder,files[i],"/3000K/STRESS")
    p=pressure.readPressureCPMD( local_file , false , 1)
    P_BLYP[i]=statistics.simpleAverage(p)
    dP_BLYP[i]=sqrt(statistics.simpleMoment(p,2))
end

V_Boates= [6.6115368589,6.7677450581,6.9401838199,7.0967427501,7.2691378687,7.426135602,7.5988814922,7.7714956737,7.86572052,7.912789266,7.9756498459,8.007058128,8.0384223944,8.0858859278,8.1013706673,8.1330423667,8.1959466238,8.3060950127,8.4943696581,8.6985238304,8.8710939961,9.0595440273,9.4521042151,9.8605002525]
P_Boates=[72.6449236765,67.3913035844,64.130432577,60.3260863203,56.8840568965,54.8913013081,52.8985522233,50.3623148819,49.275364383,48.5507242135,48.0072424605,47.6449256276,47.1014503782,48.0072424605,46.9202919617,47.6449256276,47.2826087946,46.9202919617,44.0217377873,41.6666653659,38.9492761116,36.7753621068,32.0652172641,27.7173912054]

figure()
errorbar(V,P/10,dP/10,fmt="r.-")
errorbar(Vr,Pr/10,dPr/10,fmt="b.-")
plot(V_Boates,P_Boates,"k.-")
legend(["3000K - Boates","3000K - Up","3000K - Down"])
xlabel("V (A**3/atom)")
ylabel("P (GPa)")

figure()
errorbar(V,P/10,dP/10,fmt="r.-")
errorbar(V_BLYP,P_BLYP/10,dP_BLYP/10,fmt="g.-")
plot(V_Boates,P_Boates,"k.-")
legend(["3000K - PBE - Boates","3000K - PBE","3000K - BLYP",])
xlabel("V (A**3/atom)")
ylabel("P (GPa)")

figure()
errorbar(V,P/10,dP/10,fmt="r.-")
errorbar(V_25K,P_25K/10,dP_25K/10,fmt="b.-")
errorbar(V_2K,P_2K/10,dP_2K/10,fmt="c.-")
plot(V_Boates,P_Boates,"k.-")
legend(["Boates - 3000K","3000K","2500K","2000K",])
xlabel("V (A**3/atom)")
ylabel("P (GPa)")

files=["8.82"]
local_file=string(folder,files[1],"/3000K/STRESS")
p=pressure.readPressureCPMD( local_file , false , 1)
P=Vector{Real}(Int(trunc((size(p)[1]-1000)/1000)))
Gp_size=linspace(1000,size(p)[1],1000)
for i=1:size(P)[1]
    index=Int(trunc(i*1000+1000))
    P[i]=statistics.blockAverage(p,index)
end
figure()
plot(P,"r.-")
#=========================================================#


# MSD
#=============================================================#
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
P=["8.82","9.0","9.05","9.1","9.2","9.3","9.35","9.4","9.5","9.8"]
for p in P
    file=string(folder,p,"/3000K/ENERGIES")
    MSD=CPMD.readEnergy(file)[4]
    plot(MSD)
end
legend(["8.82","9.0","9.05","9.1","9.2","9.3","9.35","9.4","9.5","9.8"])
xlabel("Time (timestep)")
ylabel("MSD(t)")

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
P=["8.82","9.0","9.4","9.8"]
for p in P
    file=string(folder,p,"/3000K/ENERGIES")
    MSD=CPMD.readEnergy(file)[4]
    plot(MSD)
end
legend(["8.82","9.0","9.4","9.8"])
xlabel("Time (timestep)")
ylabel("MSD(t)")
#=============================================================#

fileUp=string(folder,"9.4/3000K/ENERGIES")
t,e_ks,e_class,msd=CPMD.readEnergy(fileUp)
energyUp=Vector{Real}(size(e_ks)[1])
for i=1:size(e_ks)[1]
    energyUp[i] = e_ks[i]-e_class[i]
end
fileUpPress=string(folder,"9.4/3000K/STRESS")
p=pressure.readPressureCPMD( fileUpPress , false , 1)
sizep=size(p)[1]
Hu=Vector{Real}(sizep-1000)
V=9.4*9.4*9.4
for i=1001:size(p)[1]
    Hu[i-1000]=energyUp[i]+p[i]*10^8*Pascal2AU*V
end

fileDown=string(folder,"9.4r/3000K/ENERGIES")
t,e_ks,e_class,msd=CPMD.readEnergy(fileDown)
energyDown=Vector{Real}(size(e_ks)[1])
for i=1:size(e_ks)[1]
    energyDown[i] = e_ks[i]-e_class[i]
end
fileDPress=string(folder,"9.4r/3000K/STRESS")
p=pressure.readPressureCPMD( fileDPress , false , 1)
sizep=size(p)[1]
Hd=Vector{Real}(sizep)
V=9.4*9.4*9.4
for i=1:size(p)[1]
    Hd[i]=energyDown[i]+p[i]*10^8*Pascal2AU*V
end

plot(Hu)
plot(Hd)
legend(["Up","Down"])
xlabel("Time (timestep)")
ylabel("H = U + PV")
xlim([0,100000])
