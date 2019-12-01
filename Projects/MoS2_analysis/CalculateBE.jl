
# BE
#=============================================================#

# Number of Structures
Nb=[18,15,19]

# General Folder
folder="/media/moog/PNY UFD30/Data/MoS2/"

# Energy of isolated atoms
# Sulfur
file = string(folder,"S/scf.out")
S=convertRy2eV(getEnergyQE(file))
S=S[size(S)[1]]
# Mo
file = string(folder,"Mo/scf.out")
Mo =convertRy2eV(getEnergyQE(file))
Mo=Mo[size(Mo)[1]]
MoS=calBE(MoS,1,Mo,2,S)

# MoS
file = string(folder,"MoS2/relax.out")
MoS=convertRy2eV(getEnergyQE(file))
MoS=MoS[size(MoS)[1]]

# Mo2S4
Mo2S4=[]
lfold = string(folder,"MO2S4/Structures/")
for i=1:Nb[1]
  file=string(lfold,i, "/relax.out")
  conf_energy=convertRy2eV(getEnergyQE(file))
  push!(Mo2S4,conf_energy[size(conf_energy)[1]])
end
Mo2S4=calBE(Mo2S4,2,Mo,4,S)

# Mo3S6
Mo3S6=[]
lfold = string(folder,"MO3S6/Data/")
for i=1:Nb[2]
  file=string(lfold,i, "/relax.out")
  conf_energy=convertRy2eV(getEnergyQE(file))
  push!(Mo3S6,conf_energy[size(conf_energy)[1]])
end
Mo3S6=calBE(Mo3S6,3,Mo,6,S)

# Mo4S8
Mo4S8=[]
lfold = string(folder,"MO4S8/run1/Relax/")
for i=1:Nb[3]
  file=string(lfold,i, "/relax.out")
  conf_energy=convertRy2eV(getEnergyQE(file))
  push!(Mo4S8,conf_energy[size(conf_energy)[1]])
end
Mo4S8=calBE(Mo4S8,4,Mo,8,S)

# Monolayer
Mono=[]
file=string(folder,"Mono/relax.out")
conf_energy=convertRy2eV(getEnergyQE(file))
push!(Mono,conf_energy[size(conf_energy)[1]])
Mono=calBE(Mono,1,Mo,2,S)

#PW91
folderpw91="/media/moog/PNY UFD30/Data/MoS2-PW91/"
# Energy of isolated atoms
# Sulfur
filepw91 = string(folderpw91,"S/relax.out")
Spw91=convertRy2eV(getEnergyQE(filepw91))
Spw91=Spw91[size(Spw91)[1]]
# Mo
filepw91 = string(folderpw91,"Mo/relax.out")
Mopw91 =convertRy2eV(getEnergyQE(filepw91))
Mopw91=Mopw91[size(Mopw91)[1]]

# Mo3S6pw91

Mo3S6pw91=[]
lfold = string(folderpw91,"Mo3S6/")
for i=1:Nb[2]
  file=string(lfold,i, "/relax.out")
  conf_energy=convertRy2eV(getEnergyQE(file))
  push!(Mo3S6pw91,conf_energy[size(conf_energy)[1]])
end
Mo3S6pw91=calBE(Mo3S6pw91,3,Mopw91,6,Spw91)




# HOMO-LUMOgap
#==================================================================#

# Mo2S4
HL2=[]
lfold = string(folder,"MO2S4/Structures/")
for i=1:19
  file=string(lfold,i, "/relax.out")
  conf_FE=getFE(file)
  push!(HL2,conf_energy[size(conf_HL)[1]])
end
HL2=calHL(Mo2S4,2,Mo,4,S)


# Plotting
#==================================================================#
using PyPlot

BEtetra=[[1,2,3,3,3,3,4,4,4] [MoS,Mo2S4[2],Mo3S6[1],Mo3S6[5],Mo3S6[6],Mo3S6[7],Mo4S8[1],Mo4S8[4],Mo4S8[10]]]
BEtri=[[3,3,3,3,3,4,4,4,4,4,4,4,4,4] [Mo3S6[8],Mo3S6[9],Mo3S6[10],Mo3S6[12],Mo3S6[13],Mo4S8[2],Mo4S8[5],Mo4S8[7],Mo4S8[8],Mo4S8[9],Mo4S8[13],Mo4S8[15],Mo4S8[17], Mo4S8[19]]]
BEpyr=[[3,4,4,4,4,4] [Mo3S6[1],Mo4S8[1],Mo4S8[3],Mo4S8[6],Mo4S8[11],Mo4S8[12]]]
BEgraphene=[[2] [Mo2S4[17]]]
BEpent=[[2,2,3] [Mo2S4[16],Mo2S4[18],Mo2S4[15]]]
BEbulk=[[1,2,3,4] [MoS,Mo2S4[2],Mo3S6[1],Mo4S8[4]]]

# Cluster by cluster
size2=size(Mo2S4)[1]
size3=size(Mo3S6)[1]
size4=size(Mo4S8)[1]
x2=linspace(1,size2,size2)
x3=linspace(1,size3,size3)
x4=linspace(1,size4,size4)
figure()
x=linspace(0,20,20)
y=linspace(Mono[1],Mono[1],20)
plot(x,y, "c-",linewidth=3)
marker=13
plot(x2,Mo2S4, "b.-",markersize=marker)
plot(x3,Mo3S6, "g.-",markersize=marker)
plot(x4,Mo4S8, "r.-",markersize=marker)
axis([0,20,3.5,6.0])
xticks(linspace(1,19,19))
legend(["Monolayer","Mo2S4","Mo3S6","Mo4S8"],loc=3)
xlabel("# of the cluster")
ylabel("Binding Energy (eV)")
show()

# By growing n
figure()
x=linspace(0,5,20)
y=linspace(Mono[1],Mono[1],20)
x2=linspace(2,2,Nb[1])
x3=linspace(3,3,Nb[2])
x5=[]
x4=linspace(4,4,Nb[3])
marker=12
plot(x2,Mo2S4, linestyle="None",color="black", marker="." ,markersize=marker)
plot(x3,Mo3S6, linestyle="None",color="black", marker=".",markersize=marker)
plot(x4,Mo4S8, linestyle="None",color="black", marker=".",markersize=marker)
plot(BEbulk[:,1],BEbulk[:,2], linestyle="-",color="red", label="bulk-like")
plot(BEgraphene[:,1],BEgraphene[:,2], linestyle="None",color="brown", marker="." ,markersize=marker, label="graphene-like")
plot(BEpent[:,1],BEpent[:,2], linestyle="None",color="orange", marker="." ,markersize=marker, label="pentene-like")
plot(BEpyr[:,1],BEpyr[:,2], linestyfile=string(folder,"/MO2S4/MTD + Relax/Sim/COLVAR_clean")
plot(BEtri[:,1],BEtri[:,2], linestyle="None",color="green", marker="." ,markersize=marker,label="triangle")
plot(BEtetra[:,1],BEtetra[:,2], linestyle="None",color="red", marker="." ,markersize=marker,label="tetrahedron")
xticks(linspace(1,5,5))
axis([0,5,3.4,6])
legend(loc=4)
xlabel("Size of cluster (n)")
ylabel("Binding Energy (eV)")
show()


# Plot colvar
# =================================================================================



#-------------------------------------------------------------------------------

# Time Evolution of the CVs
#---------------------------
folder="/media/moog/PNY UFD30/Data/MoS2" # Work Folder
file=string(folder,"/MO2S4/MTD + Relax/Sim/COLVAR_clean") # Path to Colvar
t,cv,bias=getColvar(file,0,50000,6) # Getting information from COLVAR
t=t*0.00005; # Conversion from metasteps to ns/ps/fs

# Plot
figure()
k=0;l=0;
for i=1:nb_cvs
  if( i < 3 )
    if(k==0)
      lab="Mo"
      k=1
    else
      lab=""
    end
    plot(t,cv[:,i], color="blue",label=lab)
  else
    if(l==0)
      lab="S"
      l=1
    else
      lab=""
    end
    plot(t,cv[:,i], color="red",label=lab)
  end
end
legend()
ylim([0,4])
xlabel("Time (ps)")
ylabel("S(i) (arb. units)")
show()

#-------------------------------------------------------------------------------

# Time Evolution of the CVs + Bias potential
#----------------------------------------------
t,cv,bias=getColvar(file,0,100000,6) # Getting data from Colvar
t=t*0.00005; # Conversion of time
bias=bias;   # Conversion of the potential

# Plot
figure()
k=0;l=0;
for i=1:nb_cvs
  if( i < 3 )
    if(k==0)
      lab="Mo"
      k=1
    else
      lab=""
    end
    plot(t,cv[:,i], color="blue",label=lab)
  else
    if(l==0)
      lab="S"
      l=1
    else
      lab=""
    end
    plot(t,cv[:,i], color="red",label=lab)
  end
end
plot(t,bias,color="black",label="Bias (Ry)",linewidth=3)
ylim([0,5])
legend(loc=2)
xlabel("Time (ps)")
ylabel("Collective variables\n and Bias potential")
show()
#-------------------------------------------------------------------------------
