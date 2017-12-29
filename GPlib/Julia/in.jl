include("utils.jl");
include("cell.jl")

function readXYZ( filename::AbstractString )
  f = open(filename)
  atoms=Array(atom_basic,1)
  empty!(atoms)
  nb_atoms=parse(Int,readline(f))
  readline(f)
  if( nb_atoms > 0 )
    for i=1:nb_atoms
      line=split(readline(f))
      push!(atoms,atom_basic(str2rl(line[2]),str2rl(line[3]),str2rl(line[4]),line[1]))
    end
  elseif ( nb_atoms == 0 )
    print("Nothing to do here...\n")
    return []
  else
    error("Formatting issue with the file")
  end
  close(f)
  return atoms
end;

include("PDB.jl");
