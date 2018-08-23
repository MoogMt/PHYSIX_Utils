include("contactmatrix.jl")

# Thermodynamical values
Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]
cut_off=[1.6,1.75,1.8]

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

#for i=9:size(Volumes)[1]
for volume in Volumes
    for temperature in Temperatures
#volume=Volumes[1]
#volume=Volumes[i]
#temperature=Temperature[4]
co=cut_off[1]
# Current Volume and Temperature
current_volume=parse(Float64,volume)
current_temperature=Temperature[1]
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",current_volume,"/",temperature,"K/")
file_in=string(folder,"TRAJEC_wrapped.xyz")

# Time values

unit=0.0005*5

# Getting atoms
# for stride in strides
stride=1
atoms = filexyz.read(file_in,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

# Loop on steps to get sizes
file=open(string(folder,"atom_mol_",co,"stride",stride,".dat"),"w")
file2=open(string(folder,"largest_",co,"stride",stride,".dat"),"w")
sizemax=[]
sizes=[]
for step=1:nb_steps
    print("Building molecules: ", step/nb_steps*100,"%\n")
    #matrix = contact_matrix.buildMatrix( atoms[step] , cell, co )
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        for j=i+1:nb_atoms
            if cell_mod.distance(atoms[step],cell,i,j) < co
                matrix[i,j]=1
                matrix[j,i]=1
            end
        end
    end
    #nb_mol, mol_index = graph_mod.groupsFromMatrix(matrix,nb_atoms)
    nb_mol=0
    mol_index=zeros(nb_atoms)
    for i=1:nb_atoms
        if mol_index[i] == 0
            nb_mol += 1
            mol_index[i]=nb_mol
            mol_index = searchGroupMember(matrix,mol_index,i,nb_mol)
        end
    end
    avg_mol_size=nb_atoms/nb_mol
    # Compute sizes
    size_max=0
    for i=1:nb_mol
        size=0
        write(file,string(step," ",nb_mol," ",i," "))
        for j=1:nb_atoms
            if mol_index[j] == i
                size += 1
                write(file,string(j," "))
            end
        end
        write(file,string(size,"\n"))
        if size > size_max
            size_max=size
        end
        push!(sizes,size)
    end
    write(file2,string(step*unit*stride," ",size_max,"\n"))
end
close(file)
close(file2)

# Reading molecule file
file=open(string(folder,"atom_mol_",co,"stride",stride,".dat"))
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

check2=zeros(size(sizeVector)[1])
sum=0
for i=1:size(sizes)[1]
    check2[i] = sizes[i]*sizeVector[i]
    sum += sizes[i]*sizeVector[i]
end
check2 /= sum

file=open(string(folder,"size_proba_",co,"stride",stride,".dat"),"w")
for i=1:sizes[size(sizes)[1]]
    check_size=false
    for j=1:size(sizes)[1]
        if i == sizes[j]
            write(file,string(sizes[j]," ",check2[j],"\n"))
            chec_sizek=true
            break
        end
    end
    if check_size == false
        write(file,string(i," ",0,"\n"))
    end
end
close(file)

end
