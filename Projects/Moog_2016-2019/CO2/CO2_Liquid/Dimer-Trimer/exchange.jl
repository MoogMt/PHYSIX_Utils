# Compute the number of occurences of trimers
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

Volumes=[  9.3, 9.35, 9.375, 9.4, 9.5, 9.8, 10.0 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

#8.6, 8.82, 9.0, 9.05, 9.1, 9.15, 9.2, 9.25,

# V = 9.8
# T = 2000

nbC=32
nbO=64

cut_off = 1.75
cut_nb_= 500
max_neigh_O = 2

folder_base=string("/media/moogmt/Elements/CO2/")
folder_out = string( folder_base, "Data/Exchange/" )
if ! isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end

file_out_map_avg = string(folder_out,"avg_neighb_number_O.dat")
file_out_map_exch = string(folder_out,"exchange_number_O.dat")
handle_out_map_avg = open( file_out_map_avg, "w")
handle_out_map_exch = open( file_out_map_exch, "w")
for T in Temperatures
    file_out_map_avg_T = string(folder_out,"avg_neighb_number_O-",T,"K.dat")
    file_out_map_exch_T = string(folder_out,"exchange_number_O-",T,"K.dat")
    handle_out_map_avg_T = open( file_out_map_avg_T, "w")
    handle_out_map_exch_T = open( file_out_map_exch_T, "w")
    for V in Volumes
        folder_target = string( folder_base, V, "/", T, "K/" )
        file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )
        if ! isfile( file_traj )
            continue
        end
        traj = filexyz.readFileAtomList( file_traj )
        nb_step = size( traj )[1]
        nb_atoms = size( traj[1].names )[1]
        cell = cell_mod.Cell_param( V, V, V )
        o_neighbors = ones( nb_step, nbO, max_neigh_O )
        folder_out_local = string( folder_target, "/Data/Exchange/")
        if ! isdir( folder_out_local )
            Base.Filesystem.mkdir( folder_out_local )
        end
        file_neighbors_out = string( folder_out_local, "o_neighbors.dat")
        handle_out_neigh = open( file_neighbors_out, "w" )
        for step=1:nb_step
            print("V: ",V,"; T: ",T,"K; Progress: ",round(step/nb_step*100,digits=3),"%\n")
            for oxygen=1:nbO
                distances = zeros( nbC )
                for carbon=1:nbC
                    distances[ carbon ] = cell_mod.distance( traj[step], cell, carbon, nbC+oxygen )
                end
                index_ = sortperm( distances )
                if distances[ index_[2] ] > cut_off
                    index_[2] = -1
                end
                o_neighbors[ step, oxygen, : ] = index_[1:2]
                Base.write( file_neighbors_out, string( step, " ", oxygen, " ", index_[1], " ", index_[2], "\n" ) )
            end
        end
        close( handle_out_neigh )
        ones_v = ones(nb_step) # Dummy vector used for counting
        avg_neighbors = 0
        std_neighbors = 0
        nb_exchange = 0
        neighbors=zeros(nbO)
        for oxygen=1:nbO
            unique_neighbors = unique(o_neighbors[:,oxygen,1])
            count_actual_bonds = 0
            for i_uniq=1:size( unique_neighbors )[1]
                if unique_neighbors[i_uniq] <= 0
                    continue
                end
                nb_ = 0
                for i=1:2
                    nb_ += sum( ones_v[ o_neighbors[:,oxygen,i] .== unique_neighbors[i_uniq] ] )
                end
                if nb_ > cut_nb_
                    count_actual_bonds += 1
                end
            end
            neighbors[ oxygen ] = count_actual_bonds
            avg_neighbors += count_actual_bonds
            std_neighbors += count_actual_bonds^2
            if count_actual_bonds > 1
                nb_exchange += 1
            end
        end

        avg_neighbors /= nbO
        std_neighbors =  sqrt( std_neighbors/nbO - avg_neighbors*avg_neighbors )
        write( handle_out_map_avg, string( V," ",T," ", round(avg_neighbors,digits=3)," ",round(std_neighbors,digits=3),"\n" ) )
        write( handle_out_map_exch, string( V," ",T," ", round(nb_exchange/nbO*100,digits=3),"\n" ) )

        write( handle_out_map_avg_T, string( V," ", round(avg_neighbors,digits=3)," ",round(std_neighbors,digits=3),"\n" ) )
        write( handle_out_map_exch_T, string( V," ", round(nb_exchange/nbO*100,digits=3),"\n" ) )
    end
    close( handle_out_map_exch_T )
    close( handle_out_map_avg_T )
end
close( handle_out_map_avg )
close( handle_out_map_exch )
