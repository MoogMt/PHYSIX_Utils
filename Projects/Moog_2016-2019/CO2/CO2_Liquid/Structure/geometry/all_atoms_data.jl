# Loading necessary libraries from LibAtomicSim
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
using geom

# Objective:
# This code computes the distribution of distances and angles around carbon
# and oxygen atoms

# Sim parameters to analyze
Volumes=[10.0,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.8,8.82,8.6]
Temperatures=[1750,2000,2500,3000]

# Max number of neighbors
nbC=32
nbO=64

max_neigh=4

cut_off = 1.75
cut_off_low = 1.6

folder_base=string("/media/mathieu/Elements/CO2/")

V=8.82
T=3000

folder_target=string( folder_base, V, "/", T, "K/" )
file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )

# if isfile( file_traj )
traj = filexyz.readFileAtomList( file_traj )
species = atom_mod.getSpecies(traj[1])
species_nb = atom_mod.getNbElementSpecies( traj[1], species )

folder_out = string( folder_target, "Data/Geometry/" )
folder_out_matrix = string( folder_target, "Data/Matrix/" )

if !isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end
if !isdir( folder_out_matrix )
    Base.Filesystem.mkdir( folder_out_matrix )
end

handle_C2_X = open( string( folder_out, "C2_X.dat" ), "w" )
handle_C3_X = open( string( folder_out, "C3_X.dat" ), "w" )
handle_C4_X = open( string( folder_out, "C4_X.dat" ), "w" )
handle_C2_Y = open( string( folder_out, "C2_Y.dat" ), "w" )
handle_C3_Y = open( string( folder_out, "C3_Y.dat" ), "w" )
handle_C4_Y = open( string( folder_out, "C4_Y.dat" ), "w" )

handle_O1_X = open( string( folder_out, "O1_X.dat" ), "w" )
handle_O2_X = open( string( folder_out, "O2_X.dat" ), "w" )
handle_O1_Y = open( string( folder_out, "O1_Y.dat" ), "w" )
handle_O2_Y = open( string( folder_out, "O2_Y.dat" ), "w" )

handle_angle_C2   = open( string( folder_out, "angleC2_X.dat" ), "w" )
handle_angle_C2_Y = open( string( folder_out, "angleC2_Y.dat" ), "w" )

handle_angle_C3   = open( string( folder_out, "angleC3_X.dat" ), "w" )
handle_angle_C3_Y = open( string( folder_out, "angleC3_Y.dat" ), "w" )

handle_angle_C4   = open( string( folder_out, "angleC4_X.dat" ), "w" )
handle_angle_C4_Y = open( string( folder_out, "angleC4_Y.dat" ), "w" )

handle_angle_O2   = open( string( folder_out, "angleO2_X.dat" ), "w" )
handle_angle_O2_Y = open( string( folder_out, "angleO2_Y.dat" ), "w" )

handle_C3_short_X = open( string( folder_out, "distance_C3_short_X.dat" ), "w" )
handle_C3_long_X  = open( string( folder_out, "distance_C3_long_X.dat"  ), "w" )

handle_C3_short_Y  = open( string( folder_out, "distance_C3_short_Y.dat"  ), "w" )
handle_C3_long_Y  = open( string( folder_out, "distance_C3_long_Y.dat"  ), "w" )

handle_C4_short_X = open( string( folder_out, "distance_C4_short_X.dat" ), "w" )
handle_C4_long_X  = open( string( folder_out, "distance_C4_long_X.dat"  ), "w" )

handle_C4_short_Y = open( string( folder_out, "distance_C4_short_Y.dat" ), "w" )
handle_C4_long_Y  = open( string( folder_out, "distance_C4_long_Y.dat"  ), "w" )

handle_base_C3_X = open( string( folder_out, "base_dist_C3_X.dat"  ), "w" )
handle_base_C3_Y = open( string( folder_out, "base_dist_C3_Y.dat"  ), "w" )

handle_base_C4_X = open( string( folder_out, "base_dist_C4_X.dat"  ), "w" )
handle_base_C4_Y = open( string( folder_out, "base_dist_C4_Y.dat"  ), "w" )


handle_matrix = open( string( folder_out_matrix, "distance_matrix.dat" ), "w" )

nb_step = size(traj)[1]
nb_atoms=sum(species_nb)

cell = cell_mod.Cell_param(V,V,V)


for step=1:nb_step
    print("Progress: ",round(step/nb_step*100,digits=3),"%\n")
    distance_matrix = contact_matrix.buildMatrix( traj[step], cell)
    for carbon=1:nbC
        distances = distance_matrix[carbon,:]
        nb_neighbors = sum(ones(size(distances)[1])[ distances .< cut_off ] ) - 1
        index = sortperm(distances)
        if nb_neighbors == 1
            continue
        elseif nb_neighbors == 2
            check = true
            for neigh = 1:2
                write( handle_C2_X, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ carbon, index[2] ] # dC-O1
            b = distance_matrix[ carbon, index[3] ] # dC-O2
            c = distance_matrix[ index[2], index[3] ] # dO1-O2
            # Angle through Al-Kashi
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_C2, string( alkash_angle, "\n" ) )
            if check
                write( handle_C2_Y, string( distances[ index[ 2 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_angle_C2_Y, string( alkash_angle, "\n" ) )
            end
        elseif nb_neighbors == 3
            check = true
            for neigh = 1:3
                write( handle_C3_X, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ carbon, index[2] ] # dC-O1
            b = distance_matrix[ carbon, index[3] ] # dC-O2
            c = distance_matrix[ index[2], index[3] ] # dO1-O2
            d = distance_matrix[ carbon, index[4] ] # dC-O3
            e = distance_matrix[ index[3], index[4] ] # dO2-O3
            f = distance_matrix[ index[2], index[4] ] # dO1-O3
            # Angle through Al-Kashi
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_C2, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( a, d, f )
            write( handle_angle_C2, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( b, d, e )
            write( handle_angle_C2, string( alkash_angle, "\n" ) )
            write( handle_C3_short_X, string( distances[ index[ 2 ] ] ," " ) ) # neigh +1 because we ignore the 0
            if index[2] < nbC
                write( handle_C3_short_X, string( 0 ,"\n" ) ) # neigh +1 because we ignore the 0
            else
                write( handle_C3_short_X, string( 1 ,"\n" ) ) # neigh +1 because we ignore the 0
            end
            write( handle_C3_long_X, string( distances[ index[ 3 ] ] ,"\n" ) ) # neigh +1 because we ignore the 0
            write( handle_C3_long_X, string( distances[ index[ 4 ] ] ,"\n" ) ) # neigh +1 because we ignore the 0
            center_C = traj[step].positions[carbon,:]
            neigh_positions = zeros(Real,3,3)
            for neighb=1:3
                neigh_positions[ neighb, : ] = traj[step].positions[ index[ neighb+1 ] , : ]
                for i=1:3
                    dx = center_C[ i ] - neigh_positions[ neighb, i ]
                    if dx > V*0.5
                        neigh_positions[ neighb, i ] += V
                    end
                    if dx < -V*0.5
                        neigh_positions[ neighb, i ] -= V
                    end
                end
            end
            center_base = zeros( Real, 3 )
            for i=1:3
                center_base[i]  = ( neigh_positions[1,i] + neigh_positions[2,i] + neigh_positions[3,i] )/3
            end
            # Compute two O1-O2-O3 plan vectors
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 2, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 3, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist = dot( v_oc, normal_v )
            write( handle_base_C3_X, string( dist ,"\n") ) # Writting to disk
            if check
                write( handle_C3_Y, string( distances[ index[ 2 ] ] , "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_Y, string( distances[ index[ 3 ] ] , "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_Y, string( distances[ index[ 4 ] ] , "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_short_Y, string( distances[ index[ 2 ] ] , " " ) ) # neigh +1 because we ignore the 0
                if index[2] < nbC
                    write( handle_C3_short_Y, string( 0 ,"\n" ) ) # neigh +1 because we ignore the 0
                else
                    write( handle_C3_short_Y, string( 1 ,"\n" ) ) # neigh +1 because we ignore the 0
                end
                write( handle_C3_long_Y, string( distances[ index[ 3 ] ] , "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_long_Y, string( distances[ index[ 4 ] ] , "\n" ) ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C3_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C3_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C3_Y, string( alkash_angle, "\n" ) )
                write( handle_base_C3_Y, string( dist ,"\n") )
            end
        elseif nb_neighbors >= 4
            check = true
            for neigh = 1:4
                write( handle_C4_X, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ carbon, index[2] ] # dC-O1
            b = distance_matrix[ carbon, index[3] ] # dC-O2
            c = distance_matrix[ index[2], index[3] ] # dO1-O2
            d = distance_matrix[ carbon, index[4] ] # dC-O3
            e = distance_matrix[ index[3], index[4] ] # dO2-O3
            f = distance_matrix[ index[2], index[4] ] # dO1-O3
            g = distance_matrix[ carbon, index[5] ] # dC-O4
            h = distance_matrix[ index[3], index[5] ] # dO2-O4
            i = distance_matrix[ index[4], index[5] ] # dO3-O4
            j = distance_matrix[ index[2], index[5] ] # dO1-O4
            # Angle through Al-Kashi
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( a, d, f )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( b, d, e )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( a, g, j )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( b, g, h )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            alkash_angle = geom.angleAlKash( d, g, i )
            write( handle_angle_C4, string( alkash_angle, "\n" ) )
            # Select the target carbon positions
            center_C = traj[step].positions[carbon,:]
            # Wrap the positions of the target oxygens with regard to the positions of C
            neigh_positions = zeros(Real,4,3)
            for neighb=1:4
                neigh_positions[ neighb, : ] = traj[step].positions[ index[ neighb+1 ] , : ]
                for i=1:3
                    dx = center_C[ i ] - neigh_positions[ neighb, i ]
                    if dx > V*0.5
                        neigh_positions[ neighb, i ] += V
                    end
                    if dx < -V*0.5
                        neigh_positions[ neighb, i ] -= V
                    end
                end
            end
            # Compute two O1-O2-O3 plan vectors
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 2, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 3, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist1 = dot( v_oc, normal_v )
            write( handle_base_C4_X, string( dist1 ,"\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 3, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist2 = dot( v_oc, normal_v )
            write( handle_base_C4_X, string( dist2 ,"\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 2, : ] - neigh_positions[ 3, : ]
            v2 = neigh_positions[ 2, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            # Compute vector C-O1
            v_oc = neigh_positions[ 2, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist3 = dot( v_oc, normal_v )
            write( handle_base_C4_X, string( dist3 ,"\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 2, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist4 = dot( v_oc, normal_v )
            write( handle_base_C4_X, string( dist4 ,"\n") ) # Writting to disk
            if check
                write( handle_C4_Y, string( distances[ index[ 2 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 3 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 4 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 5 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( a, g, j )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( b, g, h )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                alkash_angle = geom.angleAlKash( d, g, i )
                write( handle_angle_C4_Y, string( alkash_angle, "\n" ) )
                write( handle_base_C4_Y, string( dist1 ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist2 ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist3 ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist4 ,"\n") ) # neigh +1 because we ignore the 0
            end
        end
    end
    for oxygen=1:nbO
        distances = distance_matrix[ nbC+oxygen,:]
        nb_neighbors = sum( ones( size( distances )[1] )[ distances .< cut_off ] ) - 1
        index = sortperm( distances )
        if nb_neighbors == 1
            check = true
            for neigh = 1:2
                write( handle_O1_X, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            if check
                write( handle_O1_Y, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
            end
        elseif nb_neighbors == 2
            check = true
            for neigh = 1:2
                write( handle_O2_X, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ nbC+oxygen, index[2] ] # dO-C1
            b = distance_matrix[ nbC+oxygen, index[3] ] # dO-C2
            c = distance_matrix[ index[2], index[3] ] # dC1-C2
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_O2, string( alkash_angle, "\n" ) )
            if check
                write( handle_O2_Y, string( distances[ index[ 2 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_O2_Y, string( distances[ index[ 3 ] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_angle_O2_Y, string( alkash_angle, "\n" ) )
            end
        else
            continue
        end
    end
    contact_matrix.writeMatrix( handle_matrix, distance_matrix )
end

close( handle_C2_X )
close( handle_C3_X )
close( handle_C4_X )

close( handle_C2_Y )
close( handle_C3_Y )
close( handle_C4_Y )

close( handle_O1_X )
close( handle_O2_X )
close( handle_O1_Y )
close( handle_O2_Y )

close( handle_angle_C2 )
close( handle_angle_C2_Y )
close( handle_angle_C3 )
close( handle_angle_C3_Y )
close( handle_angle_C4 )
close( handle_angle_C4_Y )

close( handle_angle_O2 )
close( handle_angle_O2_Y )

close( handle_C3_short_X )
close( handle_C3_long_X )

close( handle_C3_short_Y )
close( handle_C3_long_Y )

close( handle_C4_short_X )
close( handle_C4_long_X )

close( handle_C4_short_Y )
close( handle_C4_long_Y )

close( handle_base_C3_X )
close( handle_base_C3_Y )

close( handle_base_C4_X )
close( handle_base_C4_Y )

close( handle_matrix )
