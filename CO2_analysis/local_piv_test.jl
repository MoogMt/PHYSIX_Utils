GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"clustering.jl"))
include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))


# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1
min_lag=1
max_lag=5001
d_lag=5
unit=0.005


T=3000
V=8.82

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")
#folder_out=string(folder_in)

nbC=32
nbO=64
nb_atoms=nbC+nbO

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

function switchingFunction( x::T1, d::T2, n::T3, m::T4) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Int}
    return (1-(x/d)^n)/(1-(x/d)^m)
end

d=1.75
n=30
m=60

nb_type=2
step_max=100
max_neigh=5
states_sig=zeros(max_neigh*nb_type,step_max)
states_sig_sec=zeros(max_neigh*nb_type*max_neigh*nb_type,step_max)
for step=1:step_max
    bond_matrix=zeros(96,96)
    for atom1=1:nb_atoms
        for atom2=atom1+1:nb_atoms
            value=switchingFunction(cell_mod.distance(traj[step],cell,atom1,atom2),d,n,m)
            bond_matrix[atom1,atom2]=value
            bond_matrix[atom2,atom1]=value
        end
        # Sorting
        index = clustering.simpleSequence(nb_atoms)
        for carbon1=1:nbC
            for carbon2=carbon1+1:nbC
                if bond_matrix[carbon1,carbon2] < bond_matrix[carbon1,carbon3]
                    stock=bond_matrix[carbon1,carbon3]
                    bond_matrix[carbon1,carbon3] = bond_matrix[carbon1,carbon2]
                    bond_matrix[carbon1,carbon2] =stock
                    index[stock] = index[atom3]
                    index[atom3] = index[atom2]
                    index[atom2] = stock
                end
            end
        end
        for oxygen1=1:nbO
            for oxygen2=oxygen1+1:nbO
                if bond_matrix[nbC+oxygen1,nbC+oxygen2] < bond_matrix[nbC+oxygen1,nbC+oxygen3]
                    stock=bond_matrix[nbC+oxygen1,nbC+oxygen3]
                    bond_matrix[nbC+oxygen1,nbC+oxygen3] = bond_matrix[nbC+oxygen1,nbC+oxygen2]
                    bond_matrix[nbC+oxygen1,nbC+oxygen2] =stock
                    index[nbC+stock]   = index[nbC+oxygen3]
                    index[nbC+oxygen3] = index[nbC+oxygen2]
                    index[nbC+oxygen2] = stock
                end
            end
        end

    end
end
