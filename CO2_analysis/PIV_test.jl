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

phaseI, cell_I=pdb.readStep(string(folder_base,"I.pdb"))
phaseIII, cell_III=pdb.readStep(string(folder_base,"III.pdb"))

Vol_I=cell_mod.getVolume(cell_I)
Vol_III=cell_mod.getVolume(cell_III)

dist=cell_mod.distance(phaseI.positions[1,:],phaseI.positions[6,:],cell_I)

d0=5
n=4

max_=20
nb_box=200
delta=max_/nb_box
file_out=open(string(folder_base,"test_switch.dat"),"w")
for i=1:nb_box
    Base.write(file_out,string(i*delta," ",utils.switchingFunction(i*delta,d0,n),"\n"))
end
close(file_out)



nbC=864
nbO=864*2
nb_atoms=nbC+2*nbO


# positionsO_I=zeros(nbO,3)
# count_=1
# for i=1:nb_atoms
#     if phaseI.atom_names[i] == "O"
#         positionsO_I[count_,:] = phaseI.positions[i,:]
#         global count_ = count_ +  1
#     end
# end
#
# count_=1
# positionsO_III=zeros(nbO,3)
# for i=1:nb_atoms
#     if phaseIII.atom_names[i] == "O"
#         positionsO_III[count_,:] = phaseIII.positions[i,:]
#         global count_ = count_ +  1
#     end
# end
#
# count_=1
# for oxygen1=1:nbO-1
#     for oxygen2=oxygen1+1:nbO
#         piv_I[count_] = cell_mod.distance(positionsO_I[oxygen1,:],positionsO_I[oxygen2,:],cell_I)
#         piv_III[count_] = cell_mod.distance(positionsO_III[oxygen1,:],positionsO_III[oxygen2,:],cell_III)
#         global count_ = count_ + 1
#     end
# end

# ComputePIV
piv_element=Int(nbO*(nbO-1)/2)
piv_I=zeros(piv_element)
piv_III=zeros(piv_element)
count_=1
for oxygen1=1:nb_atoms
    if phaseI.atom_names[oxygen1] == "O"
        for oxygen2=oxygen1+1:nb_atoms
            if phaseI.atom_names[oxygen2] == "O"
                piv_I[count_] = utils.switchingFunction(cell_mod.distance(phaseI.positions[oxygen1,:],phaseI.positions[oxygen2,:],cell_I),d0,n)
                piv_III[count_] = utils.switchingFunction(cell_mod.distance(phaseIII.positions[oxygen1,:],phaseIII.positions[oxygen2,:],cell_III),d0,n)
                global count_ = count_ + 1
            end
        end
    end
end

# Sort
piv_I=sort(piv_I)
piv_III=sort(piv_III)

file_out=open(string(folder_base,"test_piv.dat"),"w")
for i=1:piv_element
    Base.write(file_out,string(i," ",piv_I[i]," ",piv_III[i]," ",piv_I[i]-piv_III[i],"\n"))
end
close(file_out)

d_PIV=0
for i=1:piv_element
    global d_PIV += (piv_I[i]-piv_III[i])*(piv_I[i]-piv_III[i])
end
d_PIV=sqrt(d_PIV)

print("Distance PIV = ",d_PIV,"\n")
