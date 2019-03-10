GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5


# cols (both files):
#   0:   natoms x nfi (natoms x 1, natoms x 2, ...)
#   1-3: x,y,z cartesian coords [Bohr]
#   4-6: x,y,z cartesian velocites [Bohr / thart ]
#        thart = Hartree time =  0.024189 fs
# FTRAJECTORY extra:
#   7-9: x,y,z cartesian forces [Ha / Bohr]

V=9.8
T=2500

file_in=open(string(folder_base,V,"/",T,"K/1-run/FTRAJECTORY"))
lines=readlines(file_in)
close(file_in)

nb_atoms=96
nb_steps=500

cell=cell_mod.Cell_param(V,V,V)

positions=zeros(Real,nb_steps,nb_atoms,3)
forces=zeros(Real,nb_steps,nb_atoms,3)

for i=1:nb_steps
    for j=1:nb_atoms
        line=split(lines[(i-1)*nb_atoms+j])
        for k=1:3
            positions[i,j,k] = parse( Float64, line[k+1] )
            forces[i,j,k]    = parse( Float64, line[k+7] )
        end
    end
end

max_neigh=10
distancesO=zeros(nbC*nb_steps,max_neigh)
forcesO=zeros(nbC*nb_steps,max_neigh)
for i=1:nb_steps
    for j=1:nbC
        local_distance=zeros(nbC)
        forces_local=zeros(nbC)
        for k=1:nbO
            for l=1:3
                local_distance[l] += (positions[i,j,l]-positions[i,nbC+k,l])^2
            end
            local_distance[l] = sqrt( local_distance[l] )
            forces_local[i,j,k]
        end
    end
end
