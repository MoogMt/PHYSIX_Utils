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
plot(size(sizemax)[1],sizemax,"r.")
xlabel("time (step)")
ylabel("size of largest molecule (atoms)")

# Reading molecule file
file=open(string(folder,"atoms_mol.dat"))
lines=readlines(file)
close(file)

# Keeping vectors to do size
sizes=unique(sizes)
sizeVector=zeros(size(sizes)[1])
lifetimes=[]
lifesize=[]
for index=1:size(sizes)[1]
    print(string("progress: ", index/size(sizes)[1]*100, " %\n"))
    # We get the size
    molecule_data=zeros(0,sizes[index]+2)
    for i=1:size(lines)[1]
        print(string("mini-progress: ",i/size(lines)[1]*100,"%\n"))
        line=split(lines[i])
        size_loc=parse(Int32,line[size(line)[1]])
        if sizes[index] == size_loc
            # Creating empty vector
            molecule_local=zeros(1,size_loc+2)
            # Getting step
            molecule_local[1]=parse(Int32,line[1])
            # Getting tom indexes
            for j=1:size_loc
                molecule_local[j+1]=parse(Int32,line[j+3])
            end
            # Last index is 0 (unused)
            # Adding molecule to the vector
            molecule_data=vcat(molecule_data,molecule_local)
            # Adding size occurence
            sizeVector[index] += 1
        end
    end
    # Computing lifetimes for molecule size
    nb_molecules=sizeVector[index]
    size_loc=sizes[index]
    # Loop over molecule of the same size
    for i=1:nb_molecules
        # if molecule is not used
        if molecule_data[i,size_loc+2] == 0
            # We mark it used
            molecule_data[i,size_loc+2]=1
            # Start life at 1
            life=1
            # define starting step
            start_step=molecule_data[i,1]
            # Loop over step
            for active_step=start_step+1:nb_steps
                check_step=false
                # Loop over molecules
                for k=i:nb_molecules
                    # check tool
                    check2 = false
                    # Looking at molecules that have active step
                    if molecule_data[k,size_loc+2] == active_step
                        # checking for correspondance
                        check=true
                        for l=1:size_loc
                            if molecule_data[k,l] != molecule_data[i,l]
                                check=false
                                break
                            end
                        end
                        if check
                            life += 1
                            molecule_data[k,size_loc+2] = 1
                            check_step=true
                            check2=true
                        end
                    end
                    if check2 == true
                        break
                    end
                end
                if ! check_step
                    break
                end
            end
            push!(lifesize,size_loc)
            push!(lifetimes,life)
        end
    end
end

figure(2)
plot(sizes,sizeVector,"r.")
xlabel("size (atoms)")
ylabel("Number")

check2=zeros(size(sizeVector)[1])
sum=0
for i=1:size(sizes)[1]
    check2[i] = sizes[i]*sizeVector[i]
    sum += sizes[i]*sizeVector[i]
end
check2 /= sum

figure(3)
plot(sizes,check2,"r.")
xlabel("size (atoms)")
ylabel("Frequency")

for i=1:size(lines)[1]
    line=split(lines[i])
    size_loc=parse(Int32,line[size(line)[1]])
    if 3 == size_loc
        sizeVector[index] += 1
    end
end
