GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
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

for step_sim=1:nb_steps
    bond_matrix=zeros(nb_atoms,nb_atoms)
    print("Progress: ",step_sim/nb_steps*100,"%\n")
    for atom1=1:nb_atoms
        for atom2=atom1+1:nb_atoms
            if cell_mod.distance(traj[step_sim],cell,atom1,atom2) < cut_off
                bond_matrix[atom1,atom2]=1
                bond_matrix[atom2,atom1]=1
            end
        end
    end
    coord_local=ones(nbC,max_coord,2)*(-1)
    for carbon=1:nbC
        global countC=1
        for carbon2=1:nbC
            if carbon == carbon2
                continue
            end
            # If there is a bond
            if bond_matrix[carbon,carbon2] > 0
                coord_local[carbon,countC,1] = sum(bond_matrix[carbon2,:])
                countC += 1
            end
        end
        global countO = 1
        for oxygen=1:nbO
            # If there is a bond
            if bond_matrix[carbon,oxygen+nbC] > 0
                coord_local[carbon,countO,2] = sum(bond_matrix[oxygen+nbC,:])
                countO += 1
            end
        end
    end
end
