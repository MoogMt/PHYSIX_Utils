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
            atom_mod.switchAtoms(atoms,i,j)
        end
    end
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
for i=1:count_mol
    count1=i*5+1
    if atoms.atom_names[count1] != "O1"
        for j=count1+1:5
            count2=i*5+j
            if atoms.atom_names[count2] == "O1"
                atom_mod.switchAtoms(atoms,count1,count2)
            end
        end
    end
    count1=i*5+2
    if atoms.atom_names[count1] != "C"
        for j=count1+1:5
            count2=i*5+j
            if atoms.atom_names[count2] == "C"
                atom_mod.switchAtoms(atoms,count1,count2)
            end
        end
    end
    count1=i*5+3
    if atoms.atom_names[count1] != "O2"
        for j=count1+1:5
            count2=i*5+j
            if atoms.atom_names[count2] == "O2"
                atom_mod.switchAtoms(atoms,count1,count2)
            end
        end
    end
    count1=i*5+4
    if atoms.atom_names[count1] != "M1"
        for j=count1+1:5
            count2=i*5+j
            if atoms.atom_names[count2] == "M1"
                atom_mod.switchAtoms(atoms,count1,count2)
            end
        end
    end
end
#-------------------------------------------------------------------------------

#-----------------------
# Changing atom indexes
#-------------------------------------------------------------------------------
for i=1:size(atoms.atom_names)[1]
    atoms.atom_index[i] = i
end
#-------------------------------------------------------------------------------

#-----------------
# Writting PDB
#------------------------------------------------
pdb.write(atoms,cell,"/home/moogmt/test.pdb")
#------------------------------------------------
