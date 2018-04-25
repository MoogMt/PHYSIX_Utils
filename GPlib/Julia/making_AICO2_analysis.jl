

# Loading PBD file

include("contactmatrix.jl")

folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"
atoms = filexyz.readFastFile(string(folder,"TRAJEC_wrapped.xyz"))
cell=cell_mod.Cell_param(8.82,8.82,8.82)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]


# Sensibilite to cut-off
matrix=Array{Real}(nb_atoms,nb_atoms)
mat16=Array{Real}(nb_atoms,nb_atoms)
mat17=Array{Real}(nb_atoms,nb_atoms)
mat175=Array{Real}(nb_atoms,nb_atoms)
mat18=Array{Real}(nb_atoms,nb_atoms)
for i=1:nb_steps
    matrix=contact_matrix.buildMatrix( atoms[i] , cell)

end
