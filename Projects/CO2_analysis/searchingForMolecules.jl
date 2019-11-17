# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=1
co=cut_off[1]

nbC=32
nbO=2*nbC

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/"
file=string("atom_mol_",co,"stride",stride,".dat")

if isfile( string(folder,file) )
    traj = filexyz.read(string(folder,file),stride)
    cell=cell_mod.Cell_param(V,V,V)
    nb_steps=size(traj)[1]
    nb_atoms=size(traj[1].names)[1]
    for step=1:nb_steps
        bond_matrix=zeros(nbC,nbO)
        for carbon=1:32
            for oxygen=33:96
                if cell_mod.distance(traj[step],cell,carbon,oxygen) < co
                    bond_matrix[carbon,oxygen-nbC] += 1
                end
            end
        end
    end
else
    print("file: ",file," does not exists!\n")
end
