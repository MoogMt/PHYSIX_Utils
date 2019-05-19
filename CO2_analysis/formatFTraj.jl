GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"cell.jl"))

# Folder for data
#folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

V=8.82
T=3000


file_in=open(string(folder_base,V,"/",T,"K/FTRAJECTORY2"))
folder_out=string(folder_base,V,"/",T,"K/")
lines=readlines(file_in)
close(file_in)


nb_atoms=96
nb_steps=500

cell=cell_mod.Cell_param(V,V,V)

positions=zeros(Real,nb_steps,nb_atoms,3)
forces=zeros(Real,nb_steps,nb_atoms,3)

start_step=2500

for i=start_step+1:start_step+nb_steps
    for j=1:nb_atoms
        line=split(lines[(i-1)*nb_atoms+j])
        for k=1:3
            positions[i-start_step,j,k] = cell_mod.wrap(parse( Float64, line[k+1] ),cell.length[k])
            forces[i-start_step,j,k]    = parse( Float64, line[k+7] )
        end
    end
end
