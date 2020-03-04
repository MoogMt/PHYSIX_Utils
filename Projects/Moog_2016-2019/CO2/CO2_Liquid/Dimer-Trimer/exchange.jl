# Compute the number of occurences of trimers
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

max_step=20000

Volumes=[ 8.6, 8.82, 9.0, 9.05, 9.1, 9.15, 9.2, 9.25, 9.3, 9.35, 9.375, 9.4, 9.5, 9.8, 10.0 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

V = 3000
T = 9.4

nbC=32
nbO=64

max_neigh_O = 2

folder_base=string("/media/moogmt/Elements/CO2/")
folder_out = string( folder_base, "Data/Exchange/" )
# if ! isdir( folder_out )
#     Base.Filesystem.mkdir( folder_out )
# end

# for T in Temperatures
#     for V in Volumes

folder_target = string( folder_base, V, "/", T, "K/" )

file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )

traj = filexyz.readFileAtomList( file_traj )

nb_step = size( traj )[1]
nb_atoms = size( traj[0].names )[1]

cell = cell_mod.Cell_Param( V, V, V )

o_neighbors = ones( nb_step, nbO, max_neigh_O )

for step=1:nb_step
    for oxygen=1:nbO
        distances = zeros(Int, nbC)
        for carbon=1:nbC
            distances[ carbon ] = cell_mod.distance( traj[step], cell, carbon, nbC+oxygen )
        end
        index_ = sortperm( distances )
        
    end
end


#     end
# end
