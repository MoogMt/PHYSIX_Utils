include("atoms.jl");
include("cell.jl");

module pdb

importall atom_mod
importall cell_mod

function readStep( file::AbstractString )
  #--------------
  # Reading file
  #----------------------
  file=open(file);
  lines=readlines(file);
  close(file);
  #------------------------

  #----------------------------------------------------
  # Reading informations about cell and number of atoms
  #----------------------------------------------------
  cell = cell_mod.Cell_param()
  nb_atoms=0
  for line in lines
    if split(line)[1] == "CRYST1"
      cell.a=parse(Float64,split(line)[2])
      cell.b=parse(Float64,split(line)[3])
      cell.c=parse(Float64,split(line)[4])
      cell.alpha=parse(Float64,split(line)[5])
      cell.beta=parse(Float64,split(line)[6])
      cell.gamma=parse(Float64,split(line)[7])
    elseif split(line)[1] == "ATOM"
      nb_atoms+=1
    end
  end
  #----------------------------------------------------

  #---------------------------------
  # Reading atomic informations
  #---------------------------------------------------------------------
  atoms=atom_mod.AtomMolList(nb_atoms)
  count=1
  for line in lines
    if split(line)[1] == "ATOM"
      atoms.atom_names[count]=split(line)[3]
      atoms.atom_index[count]=parse(Float64,split(line)[2])
      atoms.mol_names[count]=split(line)[4]
      atoms.mol_index[count]=parse(Float64,split(line)[5])
      atoms.positions[count,:]=[ parse(Float64,split(line)[6]), parse(Float64,split(line)[7]), parse(Float64,split(line)[8]) ]
      count+=1;
    end
  end
  #---------------------------------------------------------------------

  return atoms, cell
end

function write{ T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }(atoms::T1, cell::T2, file::T3 )

  out=open(file,"w")

  a,b,c = string(cell.a), string(cell.b), string(cell.c)
  alpha, beta, gamma = string(cell.alpha), string(cell.beta), string(cell.gamma)

  cryst1=string("CRYST1 ",a)
  cryst1=spaces(cryst1,16-length(cryst1))
  cryst1=string(cryst1,b)
  cryst1=spaces(cryst1,25-length(cryst1))
  cryst1=string(cryst1,c)
  cryst1=spaces(cryst1,34-length(cryst1))
  cryst1=string(cryst1,alpha)
  cryst1=spaces(cryst1,41-length(cryst1))
  cryst1=string(cryst1,beta)
  cryst1=spaces(cryst1,48-length(cryst1))
  cryst1=string(cryst1,gamma)
  cryst1=spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,"P 1")
  cryst1=spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  write(out,cryst1)

  write(out,string("MODEL X"))

  nb_atoms = size(atoms.atom_names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=spaces(atom,7-length(atom))
    atom=string(atom,atoms.atom_index[i])
    atom=spaces(atom,13-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=spaces(atom,23-length(atom))
    atom=string(atom,atoms.mol_names[i])
    atom=spaces(atom,27-length(atom))
    atom=string(atom,atoms.mol_index[i])
    atom=spaces(atom,31-length(atom))
    atom=string(atom,atoms.positions[i,1])
    atom=spaces(atom,39-length(atom))
    atom=string(atom,atoms.positions[i,2])
    atom=spaces(atom,47-length(atom))
    atom=string(atom,atoms.positions[i,3])
    atom=spaces(atom,55-length(atom))
    atom=string(atom,"0.00")
    atom=spaces(atom,61-length(atom))
    atom=string(atom,"0.00")
    atom=spaces(atom,77-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=string(atom,"\n")
    write(out,atom)
  end

  write(out,"END\n")

  close(out)

  return
end

end
