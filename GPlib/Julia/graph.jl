include("contactmatrix.jl")

module graph

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

function writeGraphs{ T1 <: AbstractString, T2 <: AtomList }( file::T1, Atoms_List::Vector{T2} )
    file_hand=open(file,"w")
    sizes=[]
    size_max=[]
    for step=1:nb_steps
    end
    return
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
            mol_index = graph.searchGroupMember(matrix_bonds,mol_index,i,nb_mol)
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


end
