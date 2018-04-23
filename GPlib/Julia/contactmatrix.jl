include("cell.jl");
include("xyz.jl");

module contact_matrix

using atom_mod
using cell_mod

function buildMatrix{ T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param }( atoms::T1, cell::T2 )
  nb_atoms=size(atoms.names)[1]
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    for j=1:nb_atoms
      dist=cell_mod.distance(atoms,cell,i,j)
      matrix[i,j]=dist
      matrix[j,i]=dist
    end
  end
  return matrix
end

function buildMatrix{ T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param, T3 <: Real }( atoms::T1 , cell::T2, cut_off::T3 )
  nb_atoms=size(atoms.names)[1]
  matrix=zeros(nb_atoms,nb_atoms)
  value=0
  for i=1:nb_atoms
    for j=1:nb_atoms
      dist=cell_mod.distance(atoms,cell,i,j)
      if dist <= cut_off
        matrix[i,j]=1
        matrix[j,i]=1
      else
        matrix[i,j]=0
        matrix[j,i]=0
      end
    end
  end
  return matrix
end

end
