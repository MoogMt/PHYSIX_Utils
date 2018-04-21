include("utils.jl");
include("atoms.jl");
include("cell.jl");

module contact_matrix

import atom_mod.AtomList
importall atom_mod
importall cell_mod

mutable struct ContactMatrix
    matrix::Array{Real}
    function ContactMatrix()
        new(Array{Real}(0,0))
    end
    function ContactMatrix{ T1 <: Int }( nb_atoms::T1 )
        new(zeros(nb_atoms,nb_atoms))
    end
    function ContactMatrix{ T1 <: atom_mod.AtomList }( atoms::T1 )
        nb_atoms=size(atoms)[1]
        matrix=zeros(nb_atoms,nb_atoms)
        for i=1:nb_atoms-1
            for j=i+1:nb_atoms
                dist=distance(atoms,i,j)
                matrix[i,j]=dist
                matrix[j,i]=dist
            end
        end
        new(matrix)
    end
    function ContactMatrix{ T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param }( atoms::T1 , cell::T2 )
        nb_atoms=size(atoms)[1]
        matrix=zeros(nb_atoms,nb_atoms)
        for i=1:nb_atoms-1
            for j=i+1:nb_atoms
                dist=distance(atoms,i,j)
                matrix[i,j]=dist
                matrix[j,i]=dist
            end
        end
        new(matrix)
    end
end

mutable struct ContactVector
    ContactVector::Vector{Real}
end

function BuildMatrix{ T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param }( atoms::T1, cell::T2 )
    nb_atoms=size(atoms.names)[1]
    matrix=ContactMatrix( nb_atoms )
    for i=1:nb_atoms-1
        for j=i+1:nb_atoms
            dist=distance(atoms,i,j)
            matrix.matrix[i,j]=dist
            matrix.matrix[j,i]=dist
        end
    end
    return matrix
end

end
