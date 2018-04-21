include("utils.jl");
include("atoms.jl");
include("cell.jl");

module contact_matrix

importall atom_mod
importall cell_mod

function buildMatrix{ T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param }(atoms::T1, cell::T2 )
    nb_atoms=size(atoms.names)[1]
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms-1
        for j=i+1:nb_atoms
            dist=cell_mod.distance(atoms,cell,i,j)
            matrix[i,j]=dist
            matrix[j,i]=dist
        end
    end
    return matrix
end

end
