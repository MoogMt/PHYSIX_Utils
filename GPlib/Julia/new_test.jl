include("contactmatrix.jl")

# Thermodynamical values
Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]
cut_off=[1.6,1.75,1.8]

# Current Volume and Temperature
current_volume=parse(Float64,Volumes[4])
current_temperature=Temperature[4]
folder=string("/home/moogmt/",current_volume,"/",current_temperature,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
stride=20
unit=0.0005

# Getting atoms
atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

# Loop on steps
for step=1:nb_steps
    print("Building molecules: ", step/nb_steps*100,"%\n")
    matrix = contact_matrix.buildMatrix( atoms[step] , cell, cut_off[3] )
    nb_mol, mol_index = graph_mod.groupsFromMatrix(matrix,nb_atoms)
    avg_mol_size=nb_atoms/nb_mol
    writeGraph
end

# Clearing memory
atoms=[]
