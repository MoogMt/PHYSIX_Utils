include("contactmatrix.jl")


Volumes=["8.82","9.0","9.05","9.1","9.2","9.3","9.4","9.8"]

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
folder="/home/moogmt/"
file=string(folder,"TRAJEC_wrapped.xyz")
atoms = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(8.82,8.82,8.82)
nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]
stride=5
unit=0.0005

atoms=atoms[1]

# Creating bond matrix
matrix_bonds=zeros(nb_atoms,nb_atoms)
for i=1:nb_atoms
    for j=i+1:nb_atoms
        if cell_mod.distance(atoms,cell,i,j) < 1.8
            matrix_bonds[i,j]=1
            matrix_bonds[j,i]=1
        end
    end
end
# Making molecules from bond matrix
matrix_molecules=zeros(nb_atoms,nb_atoms)
current_mol=0
nb_mol=0
for i=1:nb_atoms
    current_mol=0
    # Checking if atom is affected
    for j=1:nb_atoms
        if matrix_molecules[i,j] > 0
            current_mol=matrix_molecules[i,j]
            break
        end
    end
    # If no molecule, set
    if abs(current_mol) < 0.001
        nb_mol += 1
        current_mol= nb_mol
    end
    # Loop over each atoms
    for j=i+1:nb_atoms
        if matrix_bonds[i,j] > 0
            matrix_molecules[i,j] = current_mol
            matrix_molecules[j,i] = current_mol
        end
    end
end

f#


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

nb_box=500
dataC=makeHistogram4(timeC_hist,nb_box)
dataO=makeHistogram4(timeO_hist,nb_box)
dataC2=makeHistogram4(timeC2_hist,nb_box)
dataC3=makeHistogram4(timeC3_hist,nb_box)
dataC4=makeHistogram4(timeC4_hist,nb_box)
dataO1=makeHistogram4(timeO1_hist,nb_box)
dataO2=makeHistogram4(timeO2_hist,nb_box)
function averageVar{ T1 <: Real }( data::Array{T1} )
    avg=0
    var=0
    for i=1:size(data)[1]
        avg += data[i,1]*data[i,2]
        var += (data[i,1]^2.)*data[i,2]
    end
    return avg, var
end
datC=averageVar(dataC)
datO=averageVar(dataO)
datC2=averageVar(dataC2)
datC3=averageVar(dataC3)
datC4=averageVar(dataC4)
datO1=averageVar(dataO1)
datO2=averageVar(dataO2)
