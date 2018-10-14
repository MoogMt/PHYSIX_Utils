if ! isdefined(:cell_mod)
  include("cell.jl")
end

module filexyz

export getNbSteps, readFastFile, readStep, readEmpty, read, write

using atom_mod
using cell_mod

#------------------------------------------------------------------------------
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
function readFastFile{ T1 <: AbstractString }( file::T1 )
  #--------------
  # Reading file
  #----------------------
  file2=open(file);
  lines=readlines(file2);
  close(file2);
  #------------------------

  #------------------------
  # Basic data about files
  #-----------------------------------------
  if size(lines)[1] == 0
    print("File is empty\n")
    return sim=Vector{ atom_mod.AtomList }( 0 )
  end
  nb_atoms=parse(Int64,split(lines[1])[1])
  if size(lines)[1]/(nb_atoms+2) - trunc(size(lines)[1]/(nb_atoms+2)) > 0.00000000001
    print("ERROR: Problem inside the file! (number of atoms anounced and given do not match).")
    return sim=Vector{ atom_mod.AtomList }( 0 )
  end
  nb_steps=Int(size(lines)[1]/(nb_atoms+2))
  #------------------------------------------

  sim=Vector{ atom_mod.AtomList }( nb_steps )
  for step=1:nb_steps
      atom_list = atom_mod.AtomList( nb_atoms )
      for atom=1:nb_atoms
          line_nb=Int((step-1)*(nb_atoms+2)+atom+2)
          line_content=split( lines[line_nb] )
          atom_list.names[atom] = line_content[1]
          atom_list.index[atom] = atom
          for pos=1:3
              atom_list.positions[ atom, pos ] = parse( Float64, line_content[ pos+1 ] )
          end
      end
      sim[step]=atom_list
  end
  return sim
end
function readStep{ T1 <: IO }( file_hand::T1 )
  # Get number of atoms
  line_temp=readline(file_hand)
  nb_atoms = parse(Int,split(line_temp)[1])
  # Atoms
  atoms = AtomList(nb_atoms)
  # Reading comment line
  readline(file_hand)
  # Loop over atoms
  for i=1:nb_atoms
    line=split(readline(file_hand))
    atoms.names[i]=line[1]
    atoms.index[i]=i
    for j=1:3
      atoms.positions[i,j] = parse(Float64,line[j+1])
    end
  end
  return atoms
end
function readEmpty{ T1 <: IO }( file_hand::T1 )
  # Get number of atoms
  nb_atoms = parse(Int,split(readline(file_hand))[1])
  # Reading comment line
  readline(file_hand)
  # Loop over atoms
  for i=1:nb_atoms
    line=split(readline(file_hand))
  end
  return
end
function read{ T1 <: AbstractString }( file::T1 )
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(nb_steps)
  for i=1:nb_steps
    atoms_sim[1] = readStep(file_hand)
  end
  close(file_hand)
  return atoms_sim
end
function read{ T1 <: AbstractString, T2 <: Int }( file::T1, stride::T2 )
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(Int(trunc(nb_steps/stride)))
  j=1
  for i=1:nb_steps
    if i % stride == 0
      atoms_sim[j]=readStep(file_hand)
      j += 1
    else
      readEmpty(file_hand)
    end
  end
  return atoms_sim
end
function read{ T1 <: AbstractString, T2 <: Int , T3 <: Int }( file::T1, stride::T2 , start_step::T3 )
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(Int(trunc((nb_steps-start_step)/stride)))
  for i=1:start_step
    readEmpty(file_hand)
  end
  j=1
  for i=start_step:nb_steps-stride
    if i % stride == 0
      atoms_sim[j]=readStep(file_hand)
      j += 1
    else
      readEmpty(file_hand)
    end
  end
  return atoms_sim
end
#--------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
function write{ T1 <: IO, T2 <: atom_mod.AtomList }( file_handle::T1, atoms::T2 )
  write(file_handle,string(size(atoms.names)[1],"\n"))
  write(file_handle,string("STEP: X"))
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
  Base.write(out,string(size(atoms.names)[1],"\n"))
  Base.write(out,string("STEP: X"))
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
#---------------------------------------------------------------------------------

end
