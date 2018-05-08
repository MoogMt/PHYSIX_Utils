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

function searchGroupMember{ T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }( matrix::Array{T1}, list::Vector{T2}, index::T3 , group_nb::T4 )
    for i=1:size(matrix)[1]
        if matrix[index,i] > 0
            if list[i] == 0
                list[i]=nb_mol
                list=searchGroupMember(matrix,list,i,nb_mol)
            end
        end
    end
    return list
end

file=open(string(folder,"atoms_mol.dat"),"w")
sizes=[]
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
            mol_index=searchGroupMember(matrix_bonds,mol_index,i,nb_mol)
        end
    end

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
        size_avg += size
    end
end
close(file)

file=open(string(folder,"atoms_mol.dat"))
lines=readlines(file)
close(file)

sizes=unique(sizes)

for i=1:size(sizes)[1]
    for j=1:size(lines)[1]
        size_mol=parse(Int,lines[j])
        if size_mol == sizes[i]
            
        end
    end
end
