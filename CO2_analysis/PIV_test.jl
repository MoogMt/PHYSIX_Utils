GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using utils
using geom
using contact_matrix
using atom_mod
using cell_mod
using cube_mod
using pdb

# Folder for data
folder_base="/home/moogmt/CO2_Classic/PIV_test/"

phaseI, cell_I=pdb.readStep(string(folder_base,"I.pdb"))
phaseIII, cell_III=pdb.readStep(string(folder_base,"III.pdb"))

d0=5
n=4

nbC=864
nbO=864*2
nb_atoms=nbC+nbO

# ComputePIV
piv_element=Int(nbO*(nbO-1)/2)
piv_I=zeros(piv_element)
piv_III=zeros(piv_element)

count=1
for oxygen1=1:nb_atoms
    if phaseI.atom_names[oxygen1] != "O"
        continue
    end
    for oxygen2=oxygen1+1:nb_atoms
        if phaseI.atom_names[oxygen1] != "O"
            continue
        end
        piv_I[count] = utils.switchingFunction(cell_mod.distance(traj[step],cell,oxygen1,oxygen2),d0,n)
        piv_III[count] = utils.switchingFunction(cell_mod.distance(traj[step],cell,oxygen1,oxygen2),d0,n)
        global count = count + 1
    end
end
# Sort
piv_I=sort(piv_I)
piv_III=sort(piv_III)

d_PIV=0
for i=1:nb_structure
    d_PIV += (pivI[i]-pivIII[i])*(pivI[i]-pivIII[i])
end
d_PIV=sqrt(d_PIV)

print("Distance PIV = ",d_PIV,"\n")
