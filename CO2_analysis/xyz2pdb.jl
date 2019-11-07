GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Converting xyz trajectory to pdb trajectory

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

V=8.82
T=3000

folder=string(folder_base,V,"/",T,"K/")
folder_out=string(folder,"cluster/")

traj = filexyz.readFastFile( string(folder,"TRAJEC_wrapped.xyz") )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size( traj )[1]
nb_atoms=size( traj[1].names )[1]
nb_stride=5

file_out=open(string(folder_out,"traj.pdb"),"w")
for step=1:nb_stride:nb_steps
    truc=atom_mod.AtomList( nb_atoms )
    truc.positions=traj[step].positions
    truc.names=traj[step].names
    truc.index=traj[step].index
    pdb.write(truc,cell,file_out)
end
close(file_out)
