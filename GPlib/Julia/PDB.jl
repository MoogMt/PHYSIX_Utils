include("cell.jl");

function readPDB(file::AbstractString)
  #==========#
  # Variables
  #===========================================#
  atoms=Array(atom,1)  # Atoms
  empty!(atoms)
  #-------------------------------
  cell_name=""         # Cell name
  #--------------------------------
  cell_temp=cell([0,0,0]) # Cell structure
  #==============#
  # Reading file
  #========================================#
  open(file,"r") do f
    for line in eachline(f)
      # Getting cell information
      if ( contains(line,"CRYST1"))
        all=split(line)
        a=str2rl(all[2]); b=str2rl(all[3]); c=str2rl(all[4]);
        alpha=str2rl(all[5]); beta=str2rl(all[6]); gamma=str2rl(all[7]);
        celltype=all[8]
        cell_temp=cell(cell_param([a,b,c],[alpha,beta,gamma],celltype))
      end
      # Getting name of the cell
      if ( contains(line,"MODEL") )
        all=split(line)
        cell_name=all[2]
      end
      # Getting Atoms
      if ( contains(line,"ATOM") )
        all=split(line);
        atome_label=str2int(all[2]); atome_type=all[3];
        x=str2rl(all[6]); y=str2rl(all[7]); z=str2rl(all[8]);
        atoms=push!(atoms, atom_basic([x,y,z],atype(atome_type,atome_label)))
      end
    end
  end
  #======================#
  # Putting All Together
  #=========================================#
  # Changing cell name
  #changeCellName(cell_temp,cell_name)
  # Adding the atoms
  setAtoms(cell_temp,atoms)
  # Returning cell
  return cell_temp
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
