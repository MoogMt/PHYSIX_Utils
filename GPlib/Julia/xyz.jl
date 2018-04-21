include("utils.jl");
include("atoms.jl");
include("cell.jl");

module xyz

importall atom_mod
importall cell_mod

function getNbSteps{ T1 <: AbstractString }( file::T1)
  nb_step=0
  nb_atoms=0
  open( file ) do f
    while !eof(f)
      if nb_step == 0
        nb_atoms=parse(Float64,split(readline(f))[1])
      else
        readline(f)
      end
      nb_step+=1
    end
  end
  if nb_atoms != 0
    return Int(nb_step/(nb_atoms+2))
  else
    return 0
  end
end

function readStep{ T1 <:IO , T2 <: Int }( file_handle::T1 , nb_atoms::T2 )
  atoms = atom_mod.AtomList(nb_atoms)
  nb_atoms_check = Int(split(readline( file_handle ))[1]); readline( file_handle );
  if nb_atoms != nb_atoms_check
    print("Problem reading xyz file ! \n")
    return 0
  end
  for i=1:nb_atoms
    line=readline(file_handle)
    atoms.names[i]=split(line)[1]
    for i=1:3
      atoms.positions[i,j]=parse(Float64,split(line)[j])
    end
  end
  return atoms
end

function readFastFile{ T1 <: AbstractString }( file::T1 )
  #--------------
  # Reading file
  #----------------------
  file=open(file);
  lines=readlines(file);
  close(file);
  #------------------------

  #------------------------
  # Basic data about files
  #-----------------------------------------
  nb_atoms=parse(Int64,split(lines[1])[1])
  nb_steps=Int(size(lines)[1]/(nb_atoms+2))
  #------------------------------------------

  #-----------------------------------------------------------------------------
  atoms_boxes=Vector{ atom_mod.AtomList }(nb_steps)
  for i=1:nb_steps
    temp=atom_mod.AtomList(nb_atoms)
    for j=1:nb_atoms
      line_number=Int((i-1)*(nb_atoms+2)+j+2)
      line_content=split(lines[line_number])
      temp.names[j] = line_content[1]
      for k=1:3
        temp.positions[j,k] = parse(Float64, line_content[k+1] )
      end
      temp.index[j] = j
    end
    atoms_boxes[i] = temp
  end
  #-----------------------------------------------------------------------------

  return atoms_boxes
end

function write{ T1 <: IO, T2 <: atom_mod.AtomList }( file_handle::T1, atoms::T2 )

  for i=1:size( atoms.names )[1]
    Base.write( file_handle, string( atoms.names[i] , " "  ) )
    for j=1:3
      Base.write( file_handle, string( atoms.positions[i,j] , " " ) )
    end
    Base.write( file_handle, "\n" )
  end
end

function write{ T1 <: AbstractString, T2 <: atom_mod.AtomList }( file::T1, atoms::T2 )

  out=open(file,"w")

  for i=1:size(atoms.names)[1]
    Base.write(out, string( atoms.names[i]," " ) )
    for j=1:3
      Base.write(out, string( atoms.positions[i,j] , " ") )
    end
    Base.write(out,"\n")
  end

  close(out)
end

function write{ T1 <: AbstractString, T2 <: atom_mod.AtomList }( file::T1, atoms_blocks::Vector{T2} )
  f=open(file,"w")
  for atoms in atoms_blocks
    write(f,atoms)
  end
  close(file)
end

end
