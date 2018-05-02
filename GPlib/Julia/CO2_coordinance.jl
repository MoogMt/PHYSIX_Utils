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

function writeMatrix{ T1 <: IO , T2 <: Real }( file::T1, matrix::Array{T2} )
    nb_atoms=size(matrix)[1]
    for i=1:nb_atoms
        for j=1:nb_atoms
            write(file, string(matrix[i,j]," ") )
        end
        write(file,"\n")
    end
    return
end
stride2=2
file_matrix=open(string(folder,"coord18.dat"),"w")
coord_check=zeros(nb_atoms)
time_record=[]
time_current=zeros(nb_atoms)
for i=1:stride2:nb_steps
    coord=zeros(nb_atoms)
    write(file_matrix,string(i*stride2*stride*unit," "))
    matrix18=contact_matrix.computeMatrix(atoms[i],cell,1.8)
    for i=1:nb_atoms
        coord[i] = 0
        for j=1:nb_atoms
            coord[i] += 1
        end
        write(file_matrix,string(coord[i]," "))
    end
    avgC=0
    avgO=0
    for i=1:32
        avgC += coord[i]
    end
    avgC /= 32
    for i=33:96
        avgO += coord[i]
    end
    avgO /= 64
    write(file_matrix,string(avgC," ",avgO," "))
    write(file_matrix,"\n")
    if i > 1

    end
    coord_check=coord
end
close(file_matrix)


folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
coord=zeros(4,4,nb_steps)
file1=open(string(folder,"007_matrix_18.dat"),"w")
file2=open(string(folder,"007_matrix_175.dat"),"w")
file3=open(string(folder,"007_matrix_17.dat"),"w")
file4=open(string(folder,"007_matrix_16.dat"),"w")
write(file1,string(nb_steps," ",nb_atoms,"\n"))
write(file2,string(nb_steps," ",nb_atoms,"\n"))
write(file3,string(nb_steps," ",nb_atoms,"\n"))
write(file4,string(nb_steps," ",nb_atoms,"\n"))
matrix_bonds=zeros(nb_atoms,nb_atoms,4)
for k=1:nb_steps
    matrix_bonds=zeros(nb_atoms,nb_atoms,4)
    for j=1:96
        #-----------------------------------------------------------------------
        n_coord=zeros(4)
        # Loop over all oxygens
        for i=1:96
            #========================#
            # Distance and cut-off
            #===============================================#
            if i != j
                dist=cell_mod.distance(atoms[k],cell,j,i)
                if dist < 1.8
                    matrix_bonds[i,j,1]=1
                    matrix_bonds[j,i,1]=1
                    n_coord[1] += 1
                    if dist < 1.75
                        matrix_bonds[i,j,2]=1
                        matrix_bonds[j,i,2]=1
                        n_coord[2] += 1
                        if dist < 1.70
                            matrix_bonds[i,j,3]=1
                            matrix_bonds[j,i,3]=1
                            n_coord[3] += 1
                            if dist < 1.6
                                matrix_bonds[i,j,4]=1
                                matrix_bonds[j,i,4]=1
                                n_coord[4] += 1
                            end
                        end
                    end
                end
            end
            #===============================================#
        end
        #-----------------------------------------------------------------------
        # Computing coordinance
        #-----------------------------------------------------------------------
        for p=1:4
            if n_coord[p] == 2
                coord[p,1,k] += 1
            elseif n_coord[p] == 3
                coord[p,2,k] += 1
            elseif n_coord[p] == 4
                coord[p,3,k] += 1
            else
                coord[p,4,k] += 1
            end
        end
    end
    #-----------------------------------------------------------------------
    # Writting matrix to file
    #----------------------------------------------------------------------
    for i=1:nb_atoms
        for j=1:nb_atoms
            write(file1,string(matrix_bonds[i,j,1]," "))
            write(file1,string(matrix_bonds[i,j,2]," "))
            write(file1,string(matrix_bonds[i,j,3]," "))
            write(file1,string(matrix_bonds[i,j,4]," "))
        end
        write(file1,"\n")
        write(file2,"\n")
        write(file3,"\n")
        write(file4,"\n")
    end
    #----------------------------------------------------------------------
end
close(file1)
close(file2)
close(file3)
close(file4)

using PyPlot
figure(1)
t=linspace(0, nb_steps*stride*unit, nb_steps)
plot(t,coord[4,1,:]/32,".")
plot(t,coord[4,2,:]/32,".")
plot(t,coord[4,3,:]/32,".")
xlabel("time(ps)")
ylabel("% carbons in a given coordinance")
legend(["2C","3C","4C"])
