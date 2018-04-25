
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
coord=zeros(4,nb_steps)
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
    coord16=0; coord17=0; coord175=0; coord18=0;
    for i=1:nb_atoms
        for j=1:nb_atoms
            if matrix[i,j] < 1.6
                # mat16[i,j]=1
                # mat16[j,i]=1
                # mat17[i,j]=1
                # mat17[j,i]=1
                # mat175[i,j]=1
                # mat175[j,i]=1
                # mat18[i,j]=1
                # mat18[j,i]=1
                coord16+=1
                coord17+=1
                coord175+=1
                coord18+=1
            elseif matrix[i,j] < 1.7
                # mat17[i,j]=1
                # mat17[j,i]=1
                # mat175[i,j]=1
                # mat175[j,i]=1
                # mat18[i,j]=1
                # mat18[j,i]=1
                coord17+=1
                coord175+=1
                coord18+=1
            elseif matrix[i,j] < 1.75
                # mat175[i,j]=1
                # mat175[j,i]=1
                # mat18[i,j]=1
                # mat18[j,i]=1
                coord175+=1
                coord18+=1
            elseif matrix[i,j] < 1.8
                # mat18[i,j]=1
                # mat18[j,i]=1
                coord18+=1
            end
            Base.write(file_matrix,string(matrix[i,j]," "))
        end
        Base.write(file_matrix,"\n")
    end
    Base.write(file_matrix,"END\n")
    time[k]=k*stride*unit
    coord[1,k]=coord16/nb_atoms
    coord[2,k]=coord17/nb_atoms
    coord[3,k]=coord175/nb_atoms
    coord[4,k]=coord18/nb_atoms
    Base.write(file_coord,string(k*stride*unit," ",coord16/nb_atoms," ",coord17/nb_atoms," ",coord175/nb_atoms," ",coord18/nb_atoms,"\n"))
end
#---------------------------------------------------------

close(file_matrix)
close(file_coord)

using PyPlot

figure(1)
plot(time,coord[1])
plot(time,coord[2])
plot(time,coord[3])
plot(time,coord[4])
