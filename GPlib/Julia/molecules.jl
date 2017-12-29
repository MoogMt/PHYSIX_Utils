include("atoms.jl");
include("utils.jl");
include("list.jl");
include("atom_mol.jl");

#==================#
# Type declaration
#==============================================================================#
type molecule
  #--------------------#
  # Internal Variables
  #----------------------------------------------------------------------------#
  mol_index::Int
  atoms::Vector{atom_mol}
  #----------------------------------------------------------------------------#
  # Constructors
  #----------------------------------------------------------------------------#
  function molecule{T1 <: Int, T2 <: atom_mol}(index::T1, atoms::Vector{T2})
    new(index,atoms)
  end
  function molecule{T1 <: atom_mol, T2 <: Int}(atoms::Vector{T1}, index::T2)
    molecule(index,atoms)
  end
  function molecule{T1 <: Int, T2 <: atom_mol}(;index::T1=0, atoms::Vector{T2}=[atom_general()])
    molecule(index,atoms)
  end
  #----------------------------------------------------------------------------#
end
#==============================================================================#

#==========#
# Accesors
#==============================================================================#
# Molecules Indexes
function getMolIndex{T1 <: molecule}(molecule::T1)
  return molecule.mol_index
end
function getMolIndex{T1 <: molecule}(molecules::Vector{T1})
  indexes=Array{Int}(0)
  for molecule in molecules
    push!(indexes,getIndex(molecule))
  end
  return indexes
end
# Atoms
function getAtoms{T1 <: molecule}(molecule::T1)
  return molecule.atoms
end
function getAtoms{T1 <: molecules}(molecules::Vector{T1})
  atoms=Array{atom_mol}(0)
  for molecule in molecules
    push!(atoms,getAtoms(molecule))
  end
  return atoms
end
#==============================================================================#

#==========#
# Mutators
#==============================================================================#
# Molecules indexes
function setMolIndex{T1 <: molecule, T2 <: Int}(molecule::T1,mol_index::T2)
  molecule.mol_index = mol_index
end
function setMolIndex{T1 <: molecule, T2 <: Int}(molecules::Vector{T1},mol_indexes::Vector{T2})
  if size(molecules)[1] == size(mol_indexes)[1]
    for i=1:size(molecules)[1]
      setMolIndex(molecules[i],mol_indexes[i])
    end
  else
    error("Sizes of the molecule list and of the index list do not match!\n")
    return;
  end
end
# Atoms

#==============================================================================#

#======#
# CO2
#==============================================================================#
function createCO2{T <: atom}(atoms::Vector{T})
  molecules=Array(molecule,1)
  empty!(molecules)
  for i=1:size(atoms)[1]
    if getName(atoms[i]) == "C"
      push!(molecules,molecule([bonded()],[atoms[i]]))
      for j=1:size(atoms)[1]
        if distance(getPosition(atoms[j]),getPosition(atoms[i])) < 1.6 && i != j
          addAtom(molecules[size(molecules)[1]],atoms[j])
        end
      end
    end
  end
  return molecules
end
#==============================================================================#
