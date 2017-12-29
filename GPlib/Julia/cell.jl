include("molecules.jl");
include("utils.jl")

"""
CELL - Definition
"""

cell_element = Union{atom,molecule}

# Cell parameters
type cell_param
  # Parameters
  lengths::Vector{Real}
  angles::Vector{Real}
  ctype::AbstractString

  # Check function
  function check_celltype(cell_type::AbstractString)
    Allowed_Types=["P","A","B","C","F","I","R"]
    check=0
    for i=1:size(Allowed_Types)[1]
      if ( Allowed_Types[i] == cell_type )
        check=1
      end
    end
    if ( check == 1 )
      return true
    else
      error("Error: cell types must be P,A,B,C,F,I or R")
      return false
    end
  end
  # Inner constructorcell_elements = Union{Vector{atom},Vector{molecules}}
  function cell_param{T1 <: Real,T2 <: Real}(lengths::Vector{T1},angles::Vector{T2},ctype::AbstractString)
    if check_dim_vec(lengths,3) && check_dim_vec(angles,3) && check_celltype(ctype) && checkAngle(angles) && checkDistance(lengths)
        new(lengths,angles,ctype)
    end
  end
end

# Cell vectors
type cell_vectors
  vectors::Matrix{Real}
  function cell_vectors{T <: Real}(vectors::Matrix{T})
    if ( check_mat_dim(vectors,3,3) )
      new(vectors)
    end
  end
end

# Alternate way to declare cell_vectors
function cell_vectors{ T1 <: Real, T2 <: Real, T3 <: Real}(vectorx::Vector{T1} , vectory::Vector{T2} , vectorz::Vector{T3} )
  if ( check_dim_vec(vectorx,3) && check_dim_vec(vectory,3) && check_dim_vec(vectorz,3) )
    cell_vectors(hcat(vectorx,vectory,vectorz))
  end
end

# Definition of cell with parameters and vectors
type cell_all
  param::cell_param
  vectors::cell_vectors
end

cell_struc = Union{cell_vectors,cell_param,cell_all};


#------------------------------------------------------------------------------#
# Definition of the cell
type cell
  struc::cell_struc
  elements::Vector{cell_element}
end

#------------------------------------------------------------------------------#
# Accessors
function getType(cell::cell,atom_nb::Int)
  return cell.atoms[atom_nb].atome_type
end
# Constructors
function cell{ T1 <: cell_struc, T2 <: atom }( cell_struc::T1, atom::T2 )
  cell(cell_struc,[atom])
end
function cell{ T <: cell_struc }(cell_struc::T)
  cell(cell_struc,[])
end
function cell{ T <: Real }(abc::Vector{T})
  cell(cell_param(abc,[90,90,90],"P"),[])
end
function cell{ T1 <: Real, T2 <: atom }(abc::Vector{T1},atom::T2)
  cell(cell_param(abc,[90,90,90],"P"),[atom])
end
function cell{ T1 <: Real, T2 <: atom }(abc::Vector{T1},atoms::Vector{T2})
  cell(cell_param(abc,[90,90,90],"P"),atoms)
end

#=======#
# ATOMS
#==============================================================================#
function setAtoms{T <: atom}(cell::cell,atom::T)
  cell.atoms=atom
end
function setAtoms{T <: atom}(cell::cell,atoms::Vector{T})
  cell.elements=atoms
end
function setMolecules{T <: molecule}(cell::cell,molecules::Vector{T})
  cell.elements=molecules
end
function getElements(cell::cell)
  return cell.elements
end
function getSpecie(cell::cell,specie::AbstractString)
  atoms=[]
  for i=1:size(cell.atoms)[1]
    if ( cell.atoms[i].atom_type.name == specie )
      push!( atoms , cell.atoms[i] )
    end
  end
  return atoms
end
#==============================================================================#

# Useful cell-related functions
#-=============================================================================#
function supercell(cell::cell,nx::Int,ny::Int,nz::Int)
  n=[nx,ny,nz]
  pos=[0.,0.,0.]
  nb_atoms=size(cell.atoms)[1]
  if (typeof(cell.cell_struc)==cell_param || typeof(cell.cell_struc)==cell_all )
    for i=1:3
      cell.cell_struc.lengths[i]=cell.cell_struc.lengths[i]*n[i]
    end
    for i=1:nx
      for j=1:ny
        for k=1:nz
          for l=1:nb_atoms
            if ( i!=1 || j!=1 || k!=1 )
              pos[1] = cell.atoms[l].position[1]+i*cell.cell_struc.lengths[i]
              pos[2] = cell.atoms[l].position[2]+j*cell.cell_struc.lengths[i]
              pos[3] = cell.atoms[l].position[3]+k*cell.cell_struc.lengths[i]
              push!(cell.atoms,atom(pos[1],pos[2],pos[3],cell.atoms[l].atome_type))
            end
          end
        end
      end
    end
  else
    print("Not implemented yet...")
  end
  return cell
end
function supercell(cell::cell,n::Vector{Int})
  supercell(cell,n[1],n[2],n[3])
end
#-------------------------------------------------------------------------------
# sortAtoms
function sortAtoms(cell::cell)
  # Get the types of atoms in cell
  types=getTypes(getElements(cell))
  # Put them in different tables
  atoms=[]
  for i=1:size(types)[1]
    push!(atoms,getSpecie(cell,types[i]))
  end
  # Sort by alphabetical order
  for i=1:size(types)[1]-1
    for j=2:size(types)[1]
      if ( types[i][1] == types[j][1] )
        if( char2int(types[j][2]) > char2int(types[i][2])  )
          type_temp=types[i]
          atom_temp=atoms[i]
          types[i]=types[j]
          atom[i]=atoms[j]
          types[j]=type_temp
          atoms[i]=atom_temp
        end
      else
        if( char2int(types[j][1]) > char2int(types[i][1])  )
          type_temp=types[i]
          atom_temp=atoms[i]
          types[i]=types[j]
          atom[i]=atoms[j]
          types[j]=type_temp
          atoms[i]=atom_temp
        end
      end
    end
  end
  # Fuse The lists
  atomsnew=[]
  for i=1:size(types)[1]
    atomsnew=[atomsnew; atoms[i]]
  end
  # replace atoms
  for i=1:size(cell.atoms)[1]
    cell.atoms[i]=atomsnew[i]
  end
  return cell
end
# backInBox
#---------------------------------
# In: cell
function backInBox(cell::cell)
  for i=1:sizeof(cell.atoms)
    for j=1:3
      if ( cell.atoms.position[i] > cell.cell_struc.length[i] )
        cell.atoms.position[i]=cell.atoms.position[i]-cell.cell_struc.length[i]
      elseif ( cell.atoms.position[i] < 0 )
        cell.atoms.position[i]=cell.atoms.position[i]+cell.cell_struc.length[i]
      end
    end
  end
end
#-------------------------------------------------------------------------------
# frac2Ang
#--------------------------------------------------
# Purpose: Transforms the cells
# In: cell, the cell in which to do the conversion
# Out: cell, the modified cell
function frac2Ang(cell::cell)
  nb_atoms=size(cell.atoms)[1]
  for i=1:nb_atoms
    for j=1:3
      cell.atoms[i].position[j]=cell.atoms[i].position[j]*cell.cell_struc.lengths[j]
    end
  end
end
function translateAtoms(cell::cell,vector::Vector{Real})
  nb_atoms=size(cell.atoms)[1]
  for i=1:nb_atoms
    for j=1:3
      cell.atoms[i].position[j]=cell.atoms[i].position[j]+vector[j]
    end
  end
end
#==============================================================================#

function createCO2(cell::cell)
  molecules=Array(molecule,1)
  empty!(molecules)
  atoms=getElements(cell)
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
