
# Laio Algo
GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"xyz.jl"))
include(string(GPfolder,"cell.jl"))
include(string(GPfolder,"pdb.jl"))

# Reading PDB file
folder="/media/moogmt/Stock/CO2/Structures/Cmca/SuperFF/"
molecules, cell = pdb.readStep( string(folder,"Cmca-super.pdb") )

# C-O cut-off to build molecules
cut_off=1.4

size_molecule=3

#--------------------
# Building molecules
#------------------------------------------------------------------------------
count_mol=1
nb_atoms=size(molecules.atom_names)[1]
for atom1=1:nb_atoms
    if molecules.atom_names[atom1] == "C"
        molecules.mol_index[atom1]=count_mol
        molecules.atom_index[atom1]=(count_mol-1)*size_molecule+1
        molecules.mol_names[atom1]="COO"
        count1=1
        for atom2=1:nb_atoms
            if atom1 == atom2
                continue
            end
            dist=cell_mod.distance(molecules,cell,atom1,atom2)
            if dist < cut_off
                #Adding oxygen to molecule
                molecules.atom_names[atom2]=string("O",count1) # Changing name to O1 or O2
                molecules.mol_names[atom2]="COO"        # Changing molecule name
                molecules.mol_index[atom2]=count_mol    # Putting index of molecule
                molecules.atom_index[atom2] = molecules.atom_index[atom1] + count1
                count1 += 1
            end
        end
        global count_mol+=1
    end
end
#-------------------------------------------------------------------------------


#--------------------------------------
# Sorting atom_list by molecule index
#-------------------------------------------------------------------------------
for i=1:nb_atoms-1
    for j=i+1:nb_atoms
        if molecules.atom_index[i] > molecules.atom_index[j]
            # Storing
            atom_mod.switchAtoms(molecules,i,j)
        end
    end
end
#-------------------------------------------------------------------------------

# Moving all O atoms to 1.16 of the C
target_distance=1.163
for oxygen=1:nb_atoms
    if molecules.atom_names[oxygen] == "O1" || molecules.atom_names[oxygen] == "O2"
        for carbon=1:nb_atoms
            if molecules.atom_names[carbon] == "C"
                actual_distance=cell_mod.distance(molecules,cell,carbon,oxygen)
                if actual_distance < cut_off
                    ratio=target_distance/actual_distance
                    for i=1:3
                        molecules.positions[oxygen,i] = molecules.positions[carbon,i]+(molecules.positions[oxygen,i]-molecules.positions[carbon,i])*ratio
                    end
                end
            end
        end
    end
end


# Moving all O atoms to 1.16 of the C
for oxygen=1:nb_atoms
    if molecules.atom_names[oxygen] == "O"
        for carbon=oxygen+1:nb_atoms
            if molecules.atom_names[carbon] == "C"
                if cell_mod.distance(molecules,cell,carbon,oxygen) < 1.75
                    print("distance: ",cell_mod.distance(molecules,cell,carbon,oxygen),"\n")
                end
            end
        end
    end
end


#-----------------
# Writting PDB
#------------------------------------------------
pdb.writeStep(molecules,cell,"/home/moogmt/Cmca_ready.pdb")
#------------------------------------------------

cell_mod.wrap(molecules,cell)

pdb.writeStep(molecules,cell,"/home/moogmt/Cmca_ready_wrapped.pdb")
