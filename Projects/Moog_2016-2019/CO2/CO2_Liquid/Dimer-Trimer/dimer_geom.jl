# Compute the number of occurences of trimers
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

max_step=20000

Volumes=[ 9.375, 9.4, 9.5, 9.8, 10.0 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

min_angle=90
max_angle=180
delta_angle=1
nb_box = round(Int, ( max_angle - min_angle )/delta_angle )

folder_base=string("/media/mathieu/Elements/CO2/")
folder_out = string( folder_base, "Data/TrimersDimers/" )
if ! isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end

anglesOCO = zeros( Real, nb_box )
anglesCOC = zeros( Real, nb_box )
for T in Temperatures
    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/")
        target_traj = string( folder_in, "TRAJEC_fdb_wrapped.xyz" )
        if ! isfile( string( target_traj ) )
            continue
        end
        folder_target = string( folder_in, "/Data/Trimer/" )
        target_file = string( folder_target, "trimers_time.dat" )
        if ! isfile( target_file )
            continue
        end

        nb_lines = utils.getNbLines( target_file )

        index_C = zeros(Int, nb_lines, 2 )
        index_O = zeros(Int, nb_lines, 2 )
        steps = zeros(Int, nb_lines )

        handle_in = open( target_file )
        for line=1:nb_lines
            keywords = split( readline( handle_in ) )
            steps[ line ] = parse(Int, keywords[1] )
            for carbon=1:2
                index_C[ line, carbon ] = parse(Int, keywords[ carbon + 1 ] )
            end
            for oxygen=1:2
                index_O[ line, oxygen ] = parse(Int, keywords[ oxygen + 3 ] )
            end
        end
        close( handle_in )

        cell = cell_mod.Cell_param( V, V, V )

        traj=filexyz.readFileAtomList( target_traj )

        # Angle through Al-Kashi
        for step_=1:nb_lines
            # Angle O-C-O
            distC1O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[1], index_O[1] )
            distC1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[1], index_O[2] )
            distO1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[1], index_O[2] )
            angle_C1 = geom.angleAlKash( distC1O1, distC1O2, distO1O2)
            angleOCO[ round( Int, ( angle_C1 - min_angle )/delta_angle ) + 1 ] += 1
            distC2O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[2], index_O[1] )
            distC2O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[2], index_O[2] )
            distO1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[1], index_O[2] )
            angle_C2 = geom.angleAlKash( distC2O1, distC2O2, distO1O2)
            angleOCO[ round( Int, ( angle_C2 - min_angle )/delta_angle ) + 1 ] += 1
            # Angle C-O-C
            distC1C2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[1], index_C[2] )
            angle_O1 = geom.angleAlKash( distC1O1, distC2O1, distC1C2 )
            angle_O2 = geom.angleAlKash( distC1O2, distC2O2, distC1C2 )
            angleCOC[ round( Int, (angle_O1 - min_angle)/delta_angle ) + 1 ] += 1
            angleCOC[ round( Int, (angle_O2 - min_angle)/delta_angle ) + 1 ] += 1
        end
    end
end

for ibox = 1:nb_box
    angleCOC[ ibox ] = angleCOC[ ibox ]/( sind( ibox*delta_angle + min_angle ) )
    angleOCO[ ibox ] = angleOCO[ ibox ]/( sind( ibox*delta_angle + min_angle ) )
end
angleCOC /= sum( angleCOC )
angleOCO /= sum( angleOCO )

file_out_COC = open( string( folder_out, "angleCOC.dat" ), "w" )
file_out_OCO = open( string( folder_out, "angleOCO.dat" ), "w" )
for ibox = 1:nb_box
    Base.write( file_out_COC, string( ibox*delta_angle + min_angle, " ", angleCOC[ ibox ], "\n" ) )
    Base.write( file_out_OCO, string( ibox*delta_angle + min_angle, " ", angleOCO[ ibox ], "\n" ) )
end
close( file_out_COC )
close( file_out_OCO )
