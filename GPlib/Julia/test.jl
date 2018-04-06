include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

using PyPlot

C1=31
O1=81
O2=66
O3=70

atoms1, cell1, ELF1 = cube_mod.readCube("/home/moogmt/CO2_AIMD/ELF/0_structure/ELF.cube")
a1=cube_mod.getClosest(atoms1.positions[C1,:],ELF1)
a2=cube_mod.getClosest(atoms1.positions[O1,:],ELF1)

da=a2-a1

# Loading PBD file
include("atoms.jl")
include("cell.jl")
include("pdb.jl")

# Reading PDB file
atoms, cell=pdb.readStep( "/media/moogmt/Stock/CO2/Structures/Cmca/Conv/Cmca.pdb" )

#--------------------
# Building molecules
#------------------------------------------------------------------------------
count_mol=1
count_atoms=1
nb_atoms=size(atoms.atom_names)[1]
for i=1:nb_atoms
    if atoms.atom_names[i] == "C"
        atoms.mol_index[i]=count_mol
        atoms.mol_names[i]="CO2"
        count_OM=1
        for j=1:nb_atoms
            if atom_mod.distance(atoms.positions[i,:],atoms.positions[j,:]) < 1.6 && i != j
                #Adding oxygen to molecule
                atoms.atom_names[j]=string("O",count_OM) # Changing name to O1 or O2
                atoms.mol_names[j]="CO2"                 # Changing molecule name
                atoms.mol_index[j]=count_mol             # Putting index of molecule
                # # Adding fictif mass
                push!(atoms.mol_names,"CO2")                        # Molecule name
                push!(atoms.mol_index,count_mol)
                push!(atoms.atom_names,string("M",count_OM))
                push!(atoms.atom_index,nb_atoms+count_atoms)
                x=(atoms.positions[i,1]+atoms.positions[j,1])/2.
                y=(atoms.positions[i,2]+atoms.positions[j,2])/2.
                z=(atoms.positions[i,3]+atoms.positions[j,3])/2.
                atoms.positions=vcat(atoms.positions,[x y z])
                count_atoms+=1
                count_OM+=1
            end
        end
        count_mol+=1
    end
end
#-------------------------------------------------------------------------------

#--------------------------------------
# Sorting atom_list by molecule index
#-------------------------------------------------------------------------------
for i=1:size(atoms.atom_names)[1]-1
    for j=i+1:size(atoms.atom_names)[1]
        if atoms.mol_index[i] > atoms.mol_index[j]
            # Storing
            a_index=atoms.atom_index[i]
            a_name=atoms.atom_names[i]
            m_index=atoms.mol_index[i]
            m_name=atoms.mol_names[i]
            positions=atoms.positions[i,:]
            # Moving 1
            atoms.atom_index[i]=atoms.atom_index[j]
            atoms.atom_names[i]=atoms.atom_names[j]
            atoms.mol_index[i]=atoms.mol_index[j]
            atoms.mol_names[i]=atoms.mol_names[j]
            atoms.positions[i,:]=atoms.positions[j,:]
            # Moving 2
            atoms.atom_index[j]=a_index
            atoms.atom_names[j]=a_name
            atoms.mol_index[j]=m_index
            atoms.mol_names[j]=m_name
            atoms.positions[j,:]=positions
        end
    end
end
#-------------------------------------------------------------------------------

pdb.write(atoms,cell,"/home/moogmt/test.pdb")

#-------------------------------------------------------------------------------
file=open("/home/moogmt/test.xyz", "w")
write(file,string(size(atoms.atom_names)[1],"\n"))
write(file,"STEP X\n")
for i=1:size(atoms.atom_names)[1]
    line=string(atoms.atom_names[i]," ",atoms.positions[i,1]," ",atoms.positions[i,2]," ",atoms.positions[i,3],"\n")
    write(file,line)
end
close(file)
#-------------------------------------------------------------------------------
