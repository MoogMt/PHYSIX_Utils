GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"xyz.jl"))
include(string(GPfolder,"pdb.jl"))

# Reading PDB file
folder="/home/moogmt/II/"
molecules, cell=pdb.readStep( string(folder,"P42nmn-super.pdb") )

dCO=1.163
mo=15.9994
mc=12.011
Mt=(2*mo+mc)
Mm=Mt/2
dMM=sqrt(2*mo*Mt/Mm^2)*dCO
dCM=dMM/2
dOM=dCO-dCM

#--------------------
# Building molecules
#------------------------------------------------------------------------------
count_mol=0
count_atoms=1
nb_atoms=size(molecules.atom_names)[1]
count_check=0
count_check_2=1
for i=1:nb_atoms
    if molecules.atom_names[i] == "C"
        global count_mol+=1
        molecules.mol_index[i]=count_mol
        molecules.mol_names[i]="CO2"
        count_OM=1
        for j=1:nb_atoms
            if molecules.atom_names[j] != "O"
                continue
            end
            if cell_mod.distance(molecules,cell,i,j) < 1.6 && i != j
                molecules.mol_index[j]=count_mol
                dist=cell_mod.distance(molecules,cell,i,j)
                #Adding oxygen to molecule
                molecules.atom_names[j]=string("O",count_OM) # Changing name to O1 or O2
                molecules.mol_names[j]="CO2"                 # Changing molecule name
                molecules.mol_index[j]=count_mol
                    # Putting index of molecule
                # # Adding fictif mass
                push!(molecules.mol_names,"CO2")                        # Molecule name
                push!(molecules.mol_index,count_mol)
                push!(molecules.atom_names,string("M",count_OM))
                if count_OM > 2
                    print("C-",i,"\n")
                end
                push!(molecules.atom_index,nb_atoms+count_atoms)
                pos=zeros(Real,1,3)
                for k=1:3
                     pos[k]=molecules.positions[i,k]+(molecules.positions[i,k]-molecules.positions[j,k])/dist*dCM
                end
                molecules.positions=vcat(molecules.positions,pos)
                global count_atoms+=1
                count_OM+=1
            end
        end

    end
end
#-------------------------------------------------------------------------------

#--------------------------------------
# Sorting atom_list by molecule index
#-------------------------------------------------------------------------------
for i=1:size(molecules.atom_names)[1]-1
    for j=i+1:size(molecules.atom_names)[1]
        if molecules.mol_index[i] > molecules.mol_index[j]
            # Storing
            atom_mod.switchAtoms(molecules,i,j)
        end
    end
end
#-------------------------------------------------------------------------------

for i=1:size(molecules.atom_names)[1]
    if molecules.atom_names[i] == "M3"
        print(i,"\n")
    end
end

#---------------------------------
# Sorting atoms to match topology
#-------------------------------------------------------------------------------
for i=0:count_mol-2
    count1=5*i+1
    if molecules.atom_names[count1] != "O1"
        for j=count1+1:(i+1)*5
            count2=j
            if molecules.atom_names[count2] == "O1"
                atom_mod.switchAtoms(molecules,count1,count2)
            end
        end
    end
    count1=5*i+2
    if molecules.atom_names[count1] != "C"
        for j=count1+1:(i+1)*5
            count2=j
            if molecules.atom_names[count2] == "C"
                atom_mod.switchAtoms(molecules,count1,count2)
            end
        end
    end
    count1=5*i+3
    if molecules.atom_names[count1] != "O1"
        for j=count1+1:(i+1)*5
            count2=j
            if molecules.atom_names[count2] == "O1"
                atom_mod.switchAtoms(molecules,count1,count2)
            end
        end
    end
    count1=5*i+4
    if molecules.atom_names[count1] != "M1"
        for j=count1+1:(i+1)*5
            count2=j
            if molecules.atom_names[count2] == "M1"
                atom_mod.switchAtoms(molecules,count1,count2)
            end
        end
    end
end
#-------------------------------------------------------------------------------

#-----------------------
# Changing atom indexes
#-------------------------------------------------------------------------------
for i=1:size(molecules.atom_names)[1]
    molecules.atom_index[i] = i
end
#-------------------------------------------------------------------------------


# Computing emptyness
move=zeros(Real,3)
for i=1:size(molecules.atom_names)[1]
    for j=1:3
        if molecules.positions[i,j] < move[j]
            move[j] = molecules.positions[i,j]
        end
    end
end
for i=1:3
    move[i] = -move[i]
end
# Moving all atoms by that amount to put everything in the box
for i=1:size(molecules.atom_names)[1]
    for j=1:3
        molecules.positions[i,j] = molecules.positions[i,j]+move[j]+0.01
    end
end


#-----------------
# Writting PDB
#------------------------------------------------
pdb.write(molecules,cell,"/home/moogmt/III.pdb")
#------------------------------------------------
