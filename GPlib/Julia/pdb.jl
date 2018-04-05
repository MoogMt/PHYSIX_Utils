module pdb

include("cell.jl");


function readPDB( file::AbstractString )
  #--------------
  # Reading file
  #----------------------
  file=open(file_name);
  lines=readlines(file);
  close(file);
  #-----------------------

  #------------------------------------
  # Default values and declarations
  #--------------------------------------------------------------
  a,b,c,alpha,beta,gamma=0.,0.,0.,90.,90.,90. # Cell parameters
  #---------------------------------------------------------------

  #----------------------------------------------------
  # Reading informations
  #----------------------------------------------------
  for line in lines
    if split(line)[1] == "CRYST1"
        a=parse(Float64,split(line)[2])
        b=parse(Float64,split(line)[3])
        c=parse(Float64,split(line)[4])
        alpha=parse(Float64,split(line)[5])
        beta=parse(Float64,split(line)[6])
        gamma=parse(Float64,split(line)[7])
    elseif split(line)[1] == "ATOM"

    end
  end
  #----------------------------------------------------

  return atoms, cell
end

function writePDB(cell::cell,file::AbstractString)

  out=open(file,"w")

  a=string(cell.cell_struc.lengths[1])
  b=string(cell.cell_struc.lengths[2])
  c=string(cell.cell_struc.lengths[3])
  alpha=string(cell.cell_struc.angles[1])
  beta=string(cell.cell_struc.angles[2])
  gamma=string(cell.cell_struc.angles[3])

  c_name=cell.cell_name
  c_type=cell.cell_struc.ctype

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
  cryst1=string(cryst1,c_type," 1")
  cryst1=spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  write(out,cryst1)

  write(out,string("MODEL ",c_name))

  for i=1:size(cell.atoms)[1]
    x=string(cell.atoms[i].position[1])
    y=string(cell.atoms[i].position[2])
    z=string(cell.atoms[i].position[3])
    name=cell.atoms[i].atom_type.name
    atom="ATOM"
    atom=spaces(atom,7-length(atom))
    atom=string(atom,i)
    atom=spaces(atom,13-length(atom))
    atom=string(atom,name)
    atom=spaces(atom,23-length(atom))
    atom=string(atom,"X")
    atom=spaces(atom,27-length(atom))
    atom=string(atom,"1")
    atom=spaces(atom,31-length(atom))
    atom=string(atom,x)
    atom=spaces(atom,39-length(atom))
    atom=string(atom,y)
    atom=spaces(atom,47-length(atom))
    atom=string(atom,z)
    atom=spaces(atom,55-length(atom))
    atom=string(atom,"0.00")
    atom=spaces(atom,61-length(atom))
    atom=string(atom,"0.00")
    atom=spaces(atom,77-length(atom))
    atom=string(atom,name)
    atom=string(atom,"\n")
    write(out,atom)
  end

  write(out,"END\n")

  close(out)

  return
end

end
