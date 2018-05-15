include("contactmatrix.jl")

Volumes=["8.82","9.0","9.05","9.1","9.2","9.3","9.4","9.8"]

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
#folder="/home/moogmt/8.82/"
file=string(folder,"TRAJEC_wrapped.xyz")
atoms = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(8.82,8.82,8.82)
nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]
stride=5
unit=0.0005

function searchGroupMember{ T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }( matrix::Array{T1}, list::Vector{T2}, index::T3 , group_nb::T4 )
    for i=1:size(matrix)[1]
        if matrix[index,i] > 0
            if list[i] == 0
                list[i]=group_nb
                list=searchGroupMember(matrix,list,i,group_nb)
            end
        end
    end
    return list
end

file=open(string(folder,"atoms_mol.dat"),"w")
sizes=[]
sizemax=[]
for step=1:nb_steps
    percent=step/nb_steps
    print(string("Progres: ",percent*100," % \n"))
    # Creating bond matrix
    matrix_bonds=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        for j=i+1:nb_atoms
            if cell_mod.distance(atoms[step],cell,i,j) < 1.8
                matrix_bonds[i,j]=1
                matrix_bonds[j,i]=1
            end
        end
    end

    nb_mol=0
    mol_index=zeros(nb_atoms)
    for i=1:nb_atoms
        if mol_index[i] == 0
            nb_mol += 1
            mol_index = searchGroupMember(matrix_bonds,mol_index,i,nb_mol)
        end
    end
    size_check=0
    size_avg=0
    for i=1:nb_mol
        write(file,string(step," ",nb_mol," "))
        size=0
        write(file,string(i," "))
        for j=1:nb_atoms
            if mol_index[j] == i
                size += 1
                write(file,string(j," "))
            end
        end
        write(file,string(size," \n"))
        push!(sizes,size)
        if size > size_check
            size_check=size
        end
        size_avg += size
    end
    push!(sizemax,size_check)
end
close(file)

# Clearing memory
atoms=[]

# Plotting evolution of largest molecule as a function of time
using PyPlot
plot(1:size(sizemax)[1],sizemax,"r.")
xlabel("time (step)")
ylabel("size of largest molecule (atoms)")

# Reading molecule file
file=open(string(folder,"atoms_mol.dat"))
lines=readlines(file)
close(file)

# Keeping vectors to do size
sizes=unique(sizes)
# sorting sizes
for i=1:size(sizes)[1]
    for j=1:size(sizes)[1]
        if sizes[i] < sizes[j]
            check=sizes[i]
            sizes[i]=sizes[j]
            sizes[j]=check
        end
    end
end

sizeVector=zeros(size(sizes)[1])
file_add=string(folder,"check_data")
file_sorted=open(file_add,"w")
for index=1:size(sizes)[1]
    print(string("progress: ", index/size(sizes)[1]*100, " %\n"))
    # We get the size
    molecule_data=zeros(0,sizes[index]+2)
    for i=1:size(lines)[1]
        line=split(lines[i])
        size_loc=parse(Int32,line[size(line)[1]])
        if sizes[index] == size_loc
            # Writing step and size
            write(file_sorted,string(parse(Int32,line[1])," ",size_loc))
            # Getting tom indexes
            for j=1:size_loc
                write(file_sorted, string(" ",parse(Int32,line[j+3])) )
            end
            write(file_sorted,"\n")
            # Adding size occurence
            sizeVector[index] += 1
        end
    end
    # Computing lifetimes for molecule size
    nb_molecules=sizeVector[index]
    size_loc=sizes[index]
    # Loop over molecule of the same size
end
close(file_sorted)

# Number of occurence of each size
figure(2)
plot(sizes,sizeVector,"r.")
xlabel("size (atoms)")
ylabel("Number")

# Normalization
check2=zeros(size(sizeVector)[1])
sum=0
for i=1:size(sizes)[1]
    check2[i] = sizes[i]*sizeVector[i]
    sum += sizes[i]*sizeVector[i]
end
check2 /= sum

# Frequence at which an atom is in molecule of size X
figure(3)
plot(sizes,check2,"r.")
xlabel("size (atoms)")
ylabel("Frequency")

file=open(file_add)
lines=readlines(file)
close(file)


lifetimes=[]
lifesize=[]
offset=1
for i=1:size(sizes)[1]
    # Progress mark
    print("Progress size: ",i/size(sizes)[1]*100," %\n")
    # All
    local_size = Int(sizes[i])
    nb_molecules = Int(sizeVector[i])
    molecules=zeros(local_size,nb_molecules)
    step=zeros(nb_molecules)
    # Parsing  molecules
    for j=1:nb_molecules
        line=split(lines[Int(offset+j)])
        step[j]=parse(Int64,line[1])
        for k=1:local_size
            molecules[k,j]=parse(Float64,line[Int(k+1)])
        end
    end
    # Used vector
    used=zeros(nb_molecules)
    # Loop over molecules
    for j=1:nb_molecules
        print("Progress molecules: ",j/nb_molecules*100," %\n")
        if used[j] == 0
            # Mark molecule as used
            used[j] = 1
            # A new life is starting
            life=1
            # Determine birth step
            start_step=step[j]+1
            # Loop over steps
            for active_step=start_step:nb_steps
                check_step=false
                # Loop over molecules
                for k=1:nb_molecules
                    check_mol = false
                    # if check...
                    if active_step == step[k]
                        # Loop over atoms, check if identicals
                        for l=1:local_size
                            # Check_atoms
                            check_atoms=true
                            # Check all atoms index
                            if molecules[l,k] != molecules[l,j]
                                # if difference, mark it
                                check_atoms=false
                                # Break loop
                                break
                            end
                            # If all atoms are identical
                            if check_atoms
                                # We indicate that we found the molecule
                                check_mol = true
                                # We indicated that we have found something at that step
                                check_step = true
                                # Incrementing life
                                life += 1
                                # Indicating that atoms k is used
                                used[k] = 1
                            end
                        end
                    end
                    # If molecule is found, break and move on
                    if check_mol
                        break
                    end
                end
                # If nothing is found for a step, break
                if ! check_step
                    break
                end
            end
            push!(lifesize,local_size)
            push!(lifetimes,life)
        end
    end
    offset += nb_molecules
end
