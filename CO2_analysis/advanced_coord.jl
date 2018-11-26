GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[9.0]
Temperatures=[3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2
max_coord=5

restart=false

V=8.82
T=3000
cut_off=1.75

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_atoms=size(traj[1].names)[1]
nb_steps=size(traj)[1]

restart=true

file_coordinance=string(folder_out,"advanced_coordinance.dat")

file_out=open(file_coordinance,"w")

bond_matrix=zeros(nb_atoms,nb_atoms)

step=1

nb_type=2

bond_matrix=zeros(nb_atoms,nb_atoms)
for atom1=1:nb_atoms
    for atom2=atom1+1:nb_atoms
        if cell_mod.distance(traj[step],cell,atom1,atom2) < cut_off
            bond_matrix[atom1,atom2]=1
            bond_matrix[atom2,atom1]=1
        end
    end
end

for carbon=1:nbC
    coord_type=zeros(2)
    coord_nb=zeros(5,2)
    index_nb=zeros(5,2)
    count=0
    for atom=1:nb_atoms
        if carbon=atom
            continue
        end
        if bond_matrix[carbon,atom] > 0
            check=0
            if atom < 32
                check=1
            else
                check=2
            end
            coord_type[check] += 1
            coord_nb[check] = sum( bond_matrix[atom,:] )
            index_nb[coord_type[check],check]= atom
        end
        if sum(coord_type) > max_coord
            print("Oupsie: step-",step," carbon-",carbon," atom-",atom,"\n")
        end
    end
    write(file_coordinance,string(step," ",carbon," "))
    for type=1:nb_type
        write(file_coordinance,string(coord_type[type]," "))
        for i=1:5
            write(file_coordinance,string(coord_nb[i,type]," ",index_nb[i,type]," "))
        end
    end
    write(file_coordinance,string("\n"))
end

close(file_coordinance)
