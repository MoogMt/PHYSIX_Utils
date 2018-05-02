# Loading file
include("contactmatrix.jl")
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
file=string(folder,"TRAJEC_wrapped.xyz")
atoms = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(8.82,8.82,8.82)
nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]
stride=5
unit=0.0005

coord=zeros(4,nb_steps)
for k=1:nb_steps
    for j=1:32
        for i=1:96
            if i != j
                dist=cell_mod.distance(atoms[k],cell,j,i)
                if dist < 1.8
                    coord[1,k] += 1
                    if dist < 1.75
                        coord[2,k] += 1
                        if dist < 1.70
                            coord[3,k] += 1
                            if dist < 1.6
                                coord[4,k] += 1
                            end
                        end
                    end
                end
            end
        end
    end
    coord[:,k] /= 32
end
using PyPlot
figure(1)
t=linspace(0, nb_steps*stride*unit, nb_steps)
plot(t,coord[4,:],"-")
plot(t,coord[3,:],"-")
plot(t,coord[2,:],"-")
plot(t,coord[1,:],"-")
legend(["rc=1.6","rc=1.7","rc=1.75","rc=1.8"])
xlabel("time(ps)")
ylabel("Average Coordinance")

file=open(string(folder,"coord_avg.dat"),"w")
for i=1:size(coord[1,:])[1]
    write(file,string(t[i]))
    write(file," ")
    for j=1:4
        write(file,string(coord[j,i]))
        write(file," ")
    end
    write(file,"\n")
end
close(file)


include("contactmatrix.jl")
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
file=string(folder,"TRAJEC_wrapped.xyz")
atoms = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(8.82,8.82,8.82)
nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]
stride=5
unit=0.0005


stride2=1
file_matrix=open(string(folder,"coord18.dat"),"w")
coord_check=zeros(nb_atoms)
time_hist=[]
timeC_hist=[]
timeO_hist=[]
timeCother_hist=[]
timeC2_hist=[]
timeC3_hist=[]
timeC4_hist=[]
timeO1_hist=[]
timeO2_hist=[]
timeOother_hist=[]
time=zeros(nb_atoms)
for i=1:stride2:nb_steps
    coord=zeros(nb_atoms)
    write(file_matrix,string(i*stride2*stride*unit," "))
    # Computing BondMatrix
    matrix18=contact_matrix.computeMatrix(atoms[i],cell,1.8)
    # Computing coordinance
    #---------------------------
    for j=1:nb_atoms
        coord[j] = 0
        for k=1:nb_atoms
            if j != k
                if matrix18[j,k] > 0.0000001
                    coord[j] += 1
                end
            end
        end
        write(file_matrix,string(coord[j]," "))
    end
    #---------------------------
    # Computing averages
    #----------------------------------------------------------------------------
    avgC=0
    avgO=0
    for j=1:32
        avgC += coord[j]
    end
    avgAll=avgC
    avgC /= 32
    for j=33:96
        avgO += coord[j]
    end
    avgAll += avgO
    avgO /= 64
    avgAll /= nb_atoms
    write(file_matrix,string(avgC," ",avgO," ",avgAll,"\n"))
    #---------------------------------------------------------------------------
    if i > 1
        for j=1:nb_atoms
            if abs(coord[j] - coord_check[j]) > 0.0000001
                time[j] += 1
                time2=time[j]*stride2*stride*unit
                if j < 33
                    push!(timeC_hist,time2)
                    if coord_check[j] == 2
                        push!(timeC2_hist,time2)
                    elseif coord_check[j] == 3
                        push!(timeC3_hist,time2)
                    elseif coord_check[j] == 4
                        push!(timeC4_hist,time2)
                    else
                        push!(timeCother_hist,time2)
                    end
                else
                    push!(timeO_hist,time2)
                    if coord_check[j] == 1
                        push!(timeO1_hist,time2)
                    elseif coord_check[j] == 2
                        push!(timeO2_hist,time2)
                    else
                        push!(timeOother_hist,time2)
                    end
                end
                push!(time_hist,time2)
                time[j] = 0
            else
                time[j] += 1
            end
        end
    end
    coord_check=coord
    print("step:",i,"\n")
end
close(file_matrix)


using PyPlot
figure(1)
t=linspace(0, nb_steps*stride*unit, nb_steps)
plot(t,coord[4,1,:]/32,".")
plot(t,coord[4,2,:]/32,".")
plot(t,coord[4,3,:]/32,".")
xlabel("time(ps)")
ylabel("% carbons in a given coordinance")
legend(["2C","3C","4C"])
