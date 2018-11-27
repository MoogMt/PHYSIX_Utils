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

coord_matrix=zeros(nb_steps,nbC,2)

# file_out=open(string(folder_out,"coordinance-extended-C-",cut_off,".dat"),"w")
# for step_sim=1:nb_steps
step_sim=1
print("Progress: ",step_sim/nb_steps*100,"%\n")
# write(file_out,string(step_sim," "))
bond_matrix=zeros(nb_atoms,nb_atoms)
for atom1=1:nb_atoms
    for atom2=atom1+1:nb_atoms
        if cell_mod.distance(traj[step_sim],cell,atom1,atom2) < cut_off
            bond_matrix[atom1,atom2]=1
            bond_matrix[atom2,atom1]=1
        end
    end
end
for carbon=1:nbC
    coord_matrix[step_sim,carbon,1]=sum(bond_matrix[carbon,1:nbC])
    coord_matrix[step_sim,carbon,2]=sum(bond_matrix[carbon,nbC+1:nbC+nbO])
end
