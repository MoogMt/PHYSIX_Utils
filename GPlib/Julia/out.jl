include("cell.jl")

function writeXYZ(cell::cell,file::AbstractString)
  f = open(file, "w+")
  nb_atom=size(getElements(cell))[1]
  if ( nb_atom > 0 )
    write(f,"$nb_atom\n")
    write(f,"STEP 1\n")
    for i=1:nb_atom
      name=cell.elements[i].atom_type.name
      x=cell.elements[i].position[1]
      y=cell.elements[i].position[2]
      z=cell.elements[i].position[3]
      write(f,"$name $x $y $z\n")
    end
  else
    print("Nothing to write...")
  end
  close(f)
end

include("lammps.jl")
