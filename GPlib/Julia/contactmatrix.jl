include("cell.jl");
include("xyz.jl");

print("Loading Contact Matrix")

module contact_matrix

using atom_mod
using cell_mod

#-------------------------------------------------------------------------------
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
export buildMatrix
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function readMatrix{ T1 <: IO, T2 <: Int }( input::T1 , nb_atoms::T2)
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    line=split(readline(input))
    for j=1:nb_atoms
      matrix[i,j] = parse( Float64, line[j] )
    end
  end
  if split(readline(input))[1] == "END"
    return matrix
  else
    print("Problem while reading file")
    return
  end
end
function readMatrix{ T1 <: AbstractString }( file::T1 )
  file=open(file)
  nb_steps=Int(split(readline(file))[1])
  close(file)
  return nb_steps
end
export readMatrix
#-------------------------------------------------------------------------------

# Test purpose only, really bad at doing anything otherwise
function readMatrix{ T1<: AbstractString, T2<: Int }( file::T1, step::T2 )
  file=open(file)
  line=split(readline(file))

  nb_steps=parse(Int,line[1]);
  nb_atoms=parse(Int,line[2]);
  if nb_steps < step
    return false
  end
  for i=1:step
    for j=1:nb_atoms
      line=readline(file)
    end
  end
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    line=split(readline(file))
    for j=1:nb_atoms
      matrix[i,j]=parse(Float64,line[j])
    end
  end
  close(file)
  return matrix
end
export readMatrix

function getBonded{ T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Int , T4 <: Real }( atoms::T1, cell::T2, index::T3 , cut_off::T4 )
  nb_atoms=size(atoms.names)[1]
  for i=1:nb_atoms
    dist=cell_mod.distance(atoms,cell,index,i)
    if dist < cut_off
      print("Bond with ", i, "\n");
    end
  end
  return
end
export getBonded

function computeMatrix{ T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: Real }( atoms::T1, cell::T2 , cut_off::T3 )
  nb_atoms=size(atoms.names)[1]
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    for j=1:nb_atoms
      if i != j
        dist=cell_mod.distance(atoms,cell,i,j)
        if dist < cut_off
          matrix[i,j]=1
          matrix[j,i]=1
        end
      end
    end
  end
  return matrix
end
export computeMatrix

function writeMatrix{ T1 <: IO , T2 <: Real }( file::T1, matrix::Array{T2} )
    nb_atoms=size(matrix)[1]
    for i=1:nb_atoms
        for j=1:nb_atoms
            write(file, string(matrix[i,j]," ") )
        end
        write(file,"\n")
    end
    return
end
export writeMatrix

end
