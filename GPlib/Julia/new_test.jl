include("contactmatrix.jl")

# Thermodynamical values
Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]
cut_off=[1.6,1.75,1.8]


for volume in Volumes
# Current Volume and Temperature
current_volume=parse(Float64,volume)
current_temperature=Temperature[4]
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",current_volume,"/3000K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
strides=[1,2,5]
unit=0.0005*5

# Getting atoms
for stride in strides

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]


for co in cut_off

# Loop on steps to get sizes
file=open(string(folder,"atom_mol_",co,"stride",stride,".dat"),"w")
file2=open(string(folder,"largest_",co,"stride",stride,".dat"),"w")
sizemax=[]
for step=1:nb_steps
    print("Building molecules: ", step/nb_steps*100,"%\n")
    matrix = contact_matrix.buildMatrix( atoms[step] , cell, co )
    nb_mol, mol_index = graph_mod.groupsFromMatrix(matrix,nb_atoms)
    avg_mol_size=nb_atoms/nb_mol
    # Compute sizes
    sizes=[]
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
for i=1:size(sizes)[1]
    write(file,string(sizes[i]," ",check2[i],"\n"))
end
close(file)

file=open(file_add)
lines=readlines(file)
close(file)

lifetimes=[]
lifesize=[]
offset=0
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
            molecules[k,j]=parse(Float64,line[Int(k+2)])
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
            time=1
            # Determine birth step
            start_step=step[j]+1
            # Loop over steps
            for active_step=start_step:nb_steps
                check_step=false
                # Loop over molecules
                for k=1:nb_molecules
                    check_mol = false
                    if used[k] == 1
                        continue
                    end
                    # if check...
                    if active_step == step[k]
                        # Loop over atoms, check if identicals
                        # Check_atoms
                        check_atoms=true
                        for l=1:local_size

                            # Check all atoms index
                            if molecules[l,k] != molecules[l,j]
                                # if difference, mark it
                                check_atoms=false
                                # Break loop
                                break
                            end
                        end
                        # If all atoms are identical
                        if check_atoms
                            # We indicate that we found the molecule
                            check_mol = true
                            # We indicated that we have found something at that step
                            check_step = true
                            # Incrementing life
                            time += 1
                            # Indicating that atoms k is used
                            used[k] = 1
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
            push!(lifetimes,time)
            time=1
        end
    end
    offset += nb_molecules
end


file_life=open(string(folder,"life_all_",co,"stride",stride,".dat"),"w")
for index=1:size(sizes)[1]
    lifes=[]
    for i=1:size(lifetimes)[1]
        if lifesize[i] == sizes[index]
            push!(lifes,lifetimes[i]*unit*stride)
        end
    end
    min=0
    max=10
    N=500
    dlife=(max-min)/N
    down=min
    high=min+dlife
    center=zeros(N)
    freq=zeros(N)
    for i=1:N
        center[i]=(high+down)/2
        for j=1:size(lifes)[1]
            if lifes[j] < high && lifes[j] > down
                freq[i] += 1
            end
        end
        down += dlife
        high += dlife
    end
    # Normalization
    sum=0
    for i=1:N
    sum += freq[i]
    end
    freq /= sum
    file=open(string(folder,sizes[index],"_life.dat"),"w")
    for index2=1:N
        write(file,string(center[index2]," ",freq[index2],"\n"))
        write(file_life,string(sizes[index]," ",center[index2]," ",freq[index2],"\n"))
    end
    close(file)
    write(file_life,"\n")
end
close(file_life)
end
end
end
