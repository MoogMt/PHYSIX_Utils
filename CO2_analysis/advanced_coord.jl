GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.82]
Temperatures=[3000]
Cut_Off=[1.6]

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

file_out=open(string(folder_out,"coordinance-C-",cut_off,".dat"),"w")
for step_sim=1:nb_steps

    print("Progress: ",step_sim/nb_steps*100,"%\n")
    write(file_out,string(step_sim," "))
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
        write(file_out,string(sum(bond_matrix[carbon,1:nbC])," ",sum(bond_matrix[carbon,nbC+1:nbC+nbO])," ") )
    end
    write(file_out,string("\n"))
end
close(file_out)

global cases=ones(0,1)
global count_cases=[1]

for step_sim=1:nb_steps
    print("Progress: ",step_sim/nb_steps*100,"%\n")
    global cases
    global count_cases
    for carbon=1:nbC
        global cases
        global count_cases
        if size(cases)[1] == 0
            global cases=ones(1,2)
            cases[1,1]= coord_matrix[step_sim,carbon,1]
            cases[1,2]= coord_matrix[step_sim,carbon,2]
        else
            global new=true
            if size(cases)[2] == 1
                if cases[1] == coord_matrix[step_sim,carbon,1] && cases[2] == coord_matrix[step_sim,carbon,2]
                    global count[i] += 1
                    global new=false
                end
            else
                for i=1:size(cases)[1]
                    if cases[i,1] == coord_matrix[step_sim,carbon,1] && cases[i,2] == coord_matrix[step_sim,carbon,2]
                        global count_cases[i] += 1
                        global new=false
                    end
                end
            end
            if new
                global cases=[cases; transpose(coord_matrix[step_sim,carbon,:]) ]
                push!(count_cases,1)
            end
        end
    end
end
