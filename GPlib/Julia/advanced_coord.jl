include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.82]
Temperatures=[3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2
max_coord=4
max_coord_O=3
#
# for cut_off in Cut_Off
#     for T in Temperatures
#     for V in Volumes

V=8.82
T=3000
cut_off=1.75

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

# if ! isfile(string(folder_in,"TRAJEC_wrapped.xyz"))
#     continue
# end

folder_out=string(folder_in,"Data/")

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_atoms=size(traj[1].names)[1]
nb_steps=size(traj)[1]

file_out=open(string(folder_out,"advanced_coordinance.dat"),"w")
bond_matrix=zeros(nb_atoms,nb_atoms)
for step=1:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")

    bond_matrix=zeros(nb_atoms,nb_atoms)
    for atom1=1:nb_atoms
        for atom2=atom1+1:nb_atoms
            if cell_mod.distance(traj[step],cell,atom1,atom2) < cut_off
                bond_matrix[atom1,atom2]=1
                bond_matrix[atom2,atom1]=1
            end
        end
    end
    for carbon1=1:nbC
        write(file_out,string(step," "))
        write(file_out,string(carbon1," "))
        bond_vector=zeros(max_coord)
        count=1
        if sum(bond_matrix[carbon1,1:nbC]) > 0
            for carbon2=1:nbC
                if carbon1 == carbon2
                    continue
                end
                if bond_matrix[carbon1,carbon2] > 0
                    if count <= max_coord
                        bond_vector[count]=sum(bond_matrix[carbon2,:])
                        count+=1
                    else
                        print("Weirdness: ",step," ",carbon1," ",sum(bond_matrix[carbon1,:]),"\n")
                    end
                end
            end
        end
        bond_coord_nb=zeros(max_coord)
        for j=1:max_coord
            for i=1:max_coord
                if bond_vector[j] == i
                    bond_coord_nb[i] += 1
                end
            end
        end
        for bond=1:max_coord
            write(file_out,string(bond_coord_nb[bond]," "))
        end
        bond_vector=zeros(max_coord)
        count=1
        for oxygen=1:nbO
            if bond_matrix[carbon1,nbC+oxygen] > 0
                if count <= max_coord
                    bond_vector[count]=sum(bond_matrix[nbC+oxygen,:])
                    count+=1
                else
                    print("Weirdness: ",step," ",carbon1," ",sum(bond_matrix[carbon1,:]),"\n")
                end
            end
        end
        bond_coord_nb=zeros(max_coord_O)
        for j=1:max_coord
            for o_coord=1:max_coord_O
                if bond_vector[j] == o_coord
                    bond_coord_nb[o_coord] += 1
                end
            end
        end
        for bond=1:max_coord_O
            write(file_out,string(bond_coord_nb[bond]," "))
        end
        write(file_out,string("\n"))
    end

end
close(file_out)

#         end
#     end
# end
