module contact_matrix

import atom_mod.AtomList
importall atom_mod

mutable struct ContactMatrix
    ContactMatrix::Array{Real}
    function ContactMatrix()
        new(Array{Real}(0,0))
    end
    function ContactMatrix{ T1 <: AtomList }( atoms::T1 )
        nb_atoms=size(atoms)[1]
        matrix=zeros(nb_atoms,nb_atoms)
        for i=1:nb_atoms-1
            for j=i+1:nb_atoms
                distance=0
                for k=1:3
                    distance += (atoms.positions[i,k]-atoms.positions[j,k])^2
                end
                distance=sqrt(distance)
                matrix[i,j]=distance
                matrix[j,i]=distance
            end
        end
        new(matrix)
    end
    function ContactMatrix{ T1: AtomList, T2 <: Cell_param }( atoms::T1 , cell::T2 )
        nb_atoms=size(atoms)[1]
        matrix=zeros(nb_atoms,nb_atoms)
        for i=1:nb_atoms-1
            for j=i+1:nb_atoms
                distance=0
                for k=1:3
                    distance += (min_dist(atoms.positions[i,k]-atoms.positions[j,k], cell_param.length[k]))^2
                end
                distance=sqrt(distance)
                matrix[i,j]=distance
                matrix[j,i]=distance
            end
        end
        new(matrix)
    end
end

mutable struct ContactVector
    ContactVector::Vector{Real}
end

end
