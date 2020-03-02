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

folder_base=string("/media/mathieu/Elements/CO2/")
folder_out = string( folder_base, "Data/TrimersDimers/" )
if ! isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end
file_out_1 = open( string( folder_out, "TrimersFracMap.dat" ), "w" )
file_out_2 = open( string( folder_out, "TrimersMap.dat" ), "w" )
for T in Temperatures
    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/")
        if ! isfile( string( folder_in, "TRAJEC_fdb_wrapped.xyz" ) )
            continue
        end
        folder_target = string( folder_in, "/Data/Trimer/" )
        target_file = string( folder_target, "trimers_time.dat" )
        if ! isfile( target_file )
            continue
        end
        nb_lines = utils.getNbLines( target_file )
        hist_frac = zeros( Real, max_step )
        hist_map  = zeros( Real, max_step )
        handle_in = open( target_file )
        for line=1:nb_lines
            keyword = split( readline( handle_in ) )
            step_  = parse(Int, keyword[1] )
            hist_frac[step_]  = 1
            hist_map[step_]  += 1
        end
        close( handle_in )
        Base.write( file_out_1, string( V, " ", T, " ", sum( hist_frac )/max_step, "\n" ) )
        Base.write( file_out_2, string( V, " ", T, " ", sum( hist_map  ), "\n" ) )
    end
end
close( file_out_1 )
close( file_out_2 )
