include("contactmatrix.jl")

# Loading file
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
atoms = filexyz.readFastFile(string(folder,"TRAJEC_wrapped.xyz"))
cell=cell_mod.Cell_param(8.82,8.82,8.82)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]
stride=5
unit=0.0005

# Sensibility to cut-off
matrix=zeros(nb_atoms,nb_atoms)
time=Vector{Real}(nb_steps)
#-------------------------------------
# mat16=zeros(nb_atoms,nb_atoms)
# mat17=zeros(nb_atoms,nb_atoms)
# mat175=zeros(nb_atoms,nb_atoms)
# mat18=zeros(nb_atoms,nb_atoms)
#-------------------------------------
file_matrix=open(string(folder,"distanceMatrix.dat"),"w")
file_coord=open(string(folder,"coord_cutOff.dat"),"w")
Base.write(file_matrix,string(nb_steps," ",nb_atoms,"\n"))
#---------------------------------------------------------
for k=1:nb_steps
    matrix=contact_matrix.buildMatrix( atoms[k] , cell)
    print("timestep: ",k,"\n")
    coord=zeros(4)
    for i=1:32
        coordt=zeros(4)
        for j=1:nb_atoms
            Base.write(file_matrix,string(matrix[i,j]," "))
            if i != j
                if matrix[i,j] < 1.8
                    coordt[1] += 1
                    if matrix[i,j] < 1.75
                        coordt[2] += 1
                        if matrix[i,j] < 1.7
                            coordt[3] += 1
                            if matrix[i,j] < 1.6
                                coordt[4] += 1
                            end
                        end
                    end
                end
            end
        end
        coord += coordt
        Base.write(file_matrix,"\n")
    end
    Base.write(file_matrix,"END\n")
    time[k]=k*stride*unit
    coord /= 32
    Base.write(file_coord,string(k*stride*unit," ",coord[1]," ",coord[2]," ",coord[3]," ",coord[4],"\n"))
end
close(file_matrix)
close(file_coord)
#---------------------------------------------------------


using PyPlot

figure(1)
plot(time,coord[1,:],".-")
plot(time,coord[2,:],".-")
plot(time,coord[3,:],".-")
plot(time,coord[4,:],".-")
legend(["rc=1.80A","rc=1.75A","rc=1.70A","rc=1.60A"])
xlabel("time (ps)")
ylabel("Coordinance")


include("contactmatrix.jl")

#==============================================================================#
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
file=open(string(folder,"distanceMatrix.dat"))
line=split(readline(file))
nb_steps=parse(Int64,line[1])
nb_atoms=parse(Int64,line[2])
for i=1:nb_steps

    matrix=contact_matrix.readMatrix(file,nb_atoms)
end
close(file)
#==============================================================================#
