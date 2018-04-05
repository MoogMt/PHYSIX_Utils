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
atoms, cell=pdb.readStep( "/media/moogmt/Stock/CO2/Structures/Cmca/SuperFF/Cmca-super.pdb" )

#--------------------
# Building molecules
#------------------------------------------------------------------------------
count_mol=1
count_atoms=1
for i=1:size(atoms)[1]-1
    if atoms.atom_names[i] == "C"
        atoms.mol_index[i]=count_mol
        atoms.mol_names[i]="CO2"
        for j=i+1:size(atoms)[1]
            if distance(atoms.atom_positions[i,:],atoms.atom_positions[j,:]) < 1.6
                # Adding oxygen to molecule
                atoms.mol_names[j]="CO2"
                atoms.mol_index[j]=count_mol
                # Adding fictif mass
                push!(atoms.mol_names,"CO2")
                push!(atoms.mol_index,count_mol)
                push!(atoms.atom_names,"M1")
                push!(atoms.atom_index,size(atoms)[1]+count_atoms)
                vcat(atoms.atom_positions,[(atoms.atom_positions[i,1]-atoms.atom_positions[j,1])/2. (atoms.atom_positions[i,2]-atoms.atom_positions[j,2])/2. (atoms.atom_positions[i,3]-atoms.atom_positions[j,3])/2.])
                count_atoms+=1
            end
        end
        count_mol+=1
    end
end
#-------------------------------------------------------------------------------


# Sorting atom_list by molecule index
for i=1:size(atoms)[1]-1
    for j=1:size(atoms)[1]
        if atoms.atom_index[i] < atoms.atom_index[j]
            # Storing
            a_index=atoms.atom_index[i]
            a_name=atoms.atom_names[i]
            m_index=atoms.mol_index[i]
            m_name=atoms.mol_names[i]
            positions=atoms.atom_positions[i,:]
            # Moving 1
            atoms.atom_index[i]=atoms.atom_index[j]
            atoms.atoms_names[i]=atoms.atom_names[j]
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
