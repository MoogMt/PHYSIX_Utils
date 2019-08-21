GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using utils
using geom
#using contact_matrix
using atom_mod
using cell_mod
using cube_mod
using pdb

# Folder for data
folder_base="/home/moogmt/CO2_Classic/PIV_test/"
#folder_base="/home/moogmt/Data/PIV_test/"

phaseI,   cell_I   = pdb.readStep(string(folder_base,"I.pdb"))
phaseII,  cell_II  = pdb.readStep(string(folder_base,"II.pdb"))
phaseIII, cell_III = pdb.readStep(string(folder_base,"III.pdb"))

Vol_I   = cell_mod.getVolume( cell_I   )
Vol_II  = cell_mod.getVolume( cell_II  )
Vol_III = cell_mod.getVolume( cell_III )

d0=5
n=4

nbC=864
nbO=864*2
nb_atoms=nbC+2*nbO

# Cell matrix
cell_I_matrix   = cell_mod.params2Matrix( cell_I   )
cell_II_matrix  = cell_mod.params2Matrix( cell_II  )
cell_III_matrix = cell_mod.params2Matrix( cell_III )

# Scaling positions
positions_I_scaled   = cell_mod.getScalePosition( phaseI.positions,   cell_I_matrix   )
positions_II_scaled  = cell_mod.getScalePosition( phaseII.positions,  cell_II_matrix  )
positions_III_scaled = cell_mod.getScalePosition( phaseIII.positions, cell_III_matrix )

# ComputePIV
piv_element=Int(nbO*(nbO-1)/2)
piv_I   = zeros(piv_element)
piv_II  = zeros(piv_element)
piv_III = zeros(piv_element)
count_=1
for oxygen1=1:nb_atoms
    if phaseI.atom_names[oxygen1] == "O"
        for oxygen2=oxygen1+1:nb_atoms
            if phaseI.atom_names[oxygen2] == "O"
                piv_I[count_]    = utils.switchingFunction( (Vol_I*0.001/28)^(0.33333)*cell_mod.distanceScale( positions_I_scaled[oxygen1,:],   positions_I_scaled[oxygen2,:],   cell_I_matrix   ), d0, n )
                piv_II[count_]   = utils.switchingFunction( (Vol_I*0.001/28)^(0.33333)*cell_mod.distanceScale( positions_II_scaled[oxygen1,:],  positions_II_scaled[oxygen2,:],  cell_II_matrix  ), d0, n )
                piv_III[count_]  = utils.switchingFunction( (Vol_I*0.001/28)^(0.33333)*cell_mod.distanceScale( positions_III_scaled[oxygen1,:], positions_III_scaled[oxygen2,:], cell_III_matrix ), d0, n )
                global count_ = count_ + 1
            end
        end
    end
end

# Sort
piv_I   = sort( piv_I   )
piv_II  = sort( piv_II  )
piv_III = sort( piv_III )

file_out=open(string(folder_base,"test_piv.dat"),"w")
for i=1:piv_element
    Base.write(file_out,string(i," ",piv_I[i]," ",piv_II[i]," ",piv_III[i],"\n"))
end
close(file_out)

d_PIV_I_II   = 0
d_PIV_I_III  = 0
d_PIV_II_III = 0
for i=1:piv_element
    global d_PIV_I_II   += (piv_I[i]-piv_II[i])*(piv_I[i]-piv_II[i])
    global d_PIV_I_III  += (piv_I[i]-piv_III[i])*(piv_I[i]-piv_III[i])
    global d_PIV_II_III += (piv_II[i]-piv_III[i])*(piv_II[i]-piv_III[i])
end
d_PIV_I_II=sqrt(d_PIV_I_II)
d_PIV_I_III=sqrt(d_PIV_I_III)
d_PIV_II_III=sqrt(d_PIV_II_III)

nbC=864
nbO=864*2
nb_atoms=nbC+2*nbO

#==============================================================================#

phaseI, cell_I = pdb.readStep( string( folder_base, "I.pdb") )
phaseII, cell_II = pdb.readStep( string( folder_base, "II.pdb") )
phaseIII, cell_III = pdb.readStep( string( folder_base, "III.pdb") )

# I
#==============================================================================#
New_I=AtomMolList(nbO)
count_=1
for i=1:nb_atoms
    if phaseI.atom_names[i] == "O"
        New_I.positions[count_,:] = phaseI.positions[i,:]
        New_I.atom_index[count_] = count_
        New_I.mol_index[count_]= phaseI.mol_index[i]
        New_I.mol_names[count_] = phaseI.mol_names[i]
        New_I.atom_names[count_] = "O"
        global count_ += 1
    end
end
file_out=string(folder_base,"New_I2.pdb")
pdb.write(New_I,cell_I,file_out)
file_out=string(folder_base,"New_I_plu.pdb")
pdb.writePLUMED(New_I,cell_I,file_out)

# III
#=============================================================================#
New_III=AtomMolList(nbO)
count_=1
for i=1:nb_atoms
    if phaseI.atom_names[i] == "O"
        New_III.positions[count_,:] = phaseIII.positions[i,:]
        New_III.atom_index[count_] = count_
        New_III.mol_index[count_]= phaseIII.mol_index[i]
        New_III.mol_names[count_] = phaseIII.mol_names[i]
        New_III.atom_names[count_] = "O"
        global count_ += 1
    end
end
file_out=string(folder_base,"New_III2.pdb")
pdb.write(New_III,cell_III,file_out)
file_out=string(folder_base,"New_III_plu.pdb")
pdb.writePLUMED(New_III,cell_III,file_out)


# II
#==============================================================================#
New_II=AtomMolList(nbO)
count_=1
for i=1:nb_atoms
    if phaseI.atom_names[i] == "O"
        New_II.positions[count_,:] = phaseII.positions[i,:]
        New_II.atom_index[count_] = count_
        New_II.mol_index[count_]= phaseII.mol_index[i]
        New_II.mol_names[count_] = phaseII.mol_names[i]
        New_II.atom_names[count_] = "O"
        global count_ += 1
    end
end
file_out=string(folder_base,"New_II2.pdb")
pdb.write(New_II,cell_II,file_out)
file_out=string(folder_base,"New_II_plu.pdb")
pdb.writePLUMED(New_II,cell_II,file_out)
#==============================================================================#

print(d_PIV_I_II," ",d_PIV_I_III," ",d_PIV_II_III,"\n")

file_in=open(string(folder_base,"FRAME_TO_FRAME.MATRIX"))
lines=readlines(file_in)
close(file_in)

max_distance=parse(Float64,split(lines[1])[2])
distances_matrix=zeros(3,3)
for i=1:3
    for j=1:3
        distances_matrix[i,j]=parse(Float64,split(lines[i+1])[j])*max_distance
    end
end

d_PIV_I_II_plu   = 40.88
d_PIV_I_III_plu  = 39.47
d_PIV_II_III_plu = 23.77
