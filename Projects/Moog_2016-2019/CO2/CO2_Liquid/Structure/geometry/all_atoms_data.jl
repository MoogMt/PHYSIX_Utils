# Loading necessary libraries from LibAtomicSim
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
using geom
using utils
using conversion

using LinearAlgebra
using Statistics


# Objective:
# This code computes the distribution of distances and angles around carbon
# and oxygen atoms

function readData( file_in::T1, nb_dim::T2 ) where { T1 <: AbstractString, T2 <: Int }
    nb_lines = utils.getNbLines( file_in )
    if nb_lines == false || nb_lines == 0
        return false
    end
    handle_in = open( file_in)
    data = zeros(Real, nb_lines, nb_dim )
    for line=1:nb_lines
        keywords = split( readline( handle_in ) )
        for i=1:nb_dim
            data[ line, i ] = parse(Float64, keywords[i] )
        end
    end
    close( handle_in )
    return data
end

function makeHist( data::Array{T1,2}, nb_box_hist::T2, max_step::T3, block_size::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int }
    min_ = minimum( data[:,1] )
    max_ = maximum( data[:,1])
    delta_ = (max_ - min_)/nb_box_hist
    nb_block = round(Int,max_step/block_size)
    hist_boxes = zeros(Real, nb_box_hist+1, nb_block+1 )
    size_data=size(data)[2]
    for i=1:size(data[:,1])[1]
        hist_boxes[ round(Int,( data[i,1] - min_ )/delta_ )+1, Int( trunc( data[i,size_data]/block_size ) )+1  ] += 1
    end
    for i_block=1:nb_block+1
         sum_  = sum( hist_boxes[ :, i_block ] )
         #print("sum: ",sum_,"\n")
         for i_box=1:nb_box_hist+1
             hist_boxes[ i_box, i_block ] = hist_boxes[ i_box, i_block ]/sum_
         end
    end
    hist_avg = zeros( nb_box_hist+1 )
    hist_std = zeros( nb_box_hist+1 )
    for i_box=1:nb_box_hist+1
        hist_avg[ i_box ] = Statistics.mean( hist_boxes[ i_box, : ] )
        hist_std[ i_box ] = Statistics.std(   hist_boxes[ i_box, : ] )
    end
    return hist_avg, hist_std, delta_, min_
end

function writeHist( file_out::T1, hist_avg::Vector{T2}, hist_err::Vector{T3} , delta_::T4, min_::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real }
    handle_out = open( file_out, "w" )
    nb_box=size(hist_avg)[1]
    for i_box=1:nb_box
        Base.write( handle_out, string( i_box*delta_+min_, " ", hist_avg[ i_box ], " ", hist_err[ i_box ], "\n" ) )
    end
    close( handle_out )
    return true
end



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

handle_angle_C2_X = open( string( folder_out, "angleC2_X.dat" ), "w" )
handle_angle_C2_Y = open( string( folder_out, "angleC2_Y.dat" ), "w" )

handle_angle_C3   = open( string( folder_out, "angleC3_X.dat" ), "w" )
handle_angle_C3_Y = open( string( folder_out, "angleC3_Y.dat" ), "w" )

handle_angle_C3_shortlong_X   = open( string( folder_out, "angleC3_X_shortlong.dat" ), "w" )
handle_angle_C3_longlong_X    = open( string( folder_out, "angleC3_X_longlong.dat" ), "w" )

handle_angle_C3_shortlong_Y   = open( string( folder_out, "angleC3_Y_shortlong.dat" ), "w" )
handle_angle_C3_longlong_Y    = open( string( folder_out, "angleC3_Y_longlong.dat" ), "w" )

handle_angle_C3_23_X   = open( string( folder_out, "angleC3_23_X.dat" ), "w" )
handle_angle_C3_24_X   = open( string( folder_out, "angleC3_24_X.dat" ), "w" )
handle_angle_C3_34_X   = open( string( folder_out, "angleC3_34_X.dat" ), "w" )

handle_angle_C3_23_Y   = open( string( folder_out, "angleC3_23_Y.dat" ), "w" )
handle_angle_C3_24_Y   = open( string( folder_out, "angleC3_24_Y.dat" ), "w" )
handle_angle_C3_34_Y   = open( string( folder_out, "angleC3_34_Y.dat" ), "w" )

handle_angle_C4   = open( string( folder_out, "angleC4_X.dat" ), "w" )
handle_angle_C4_Y = open( string( folder_out, "angleC4_Y.dat" ), "w" )

handle_angle_C4_23_X = open( string( folder_out, "angleC4_23_X.dat" ), "w" )
handle_angle_C4_24_X = open( string( folder_out, "angleC4_24_X.dat" ), "w" )
handle_angle_C4_25_X = open( string( folder_out, "angleC4_25_X.dat" ), "w" )
handle_angle_C4_34_X = open( string( folder_out, "angleC4_34_X.dat" ), "w" )
handle_angle_C4_35_X = open( string( folder_out, "angleC4_35_X.dat" ), "w" )
handle_angle_C4_45_X = open( string( folder_out, "angleC4_45_X.dat" ), "w" )

handle_angle_C4_23_Y = open( string( folder_out, "angleC4_23_Y.dat" ), "w" )
handle_angle_C4_24_Y = open( string( folder_out, "angleC4_24_Y.dat" ), "w" )
handle_angle_C4_25_Y = open( string( folder_out, "angleC4_25_Y.dat" ), "w" )
handle_angle_C4_34_Y = open( string( folder_out, "angleC4_34_Y.dat" ), "w" )
handle_angle_C4_35_Y = open( string( folder_out, "angleC4_35_Y.dat" ), "w" )
handle_angle_C4_45_Y = open( string( folder_out, "angleC4_45_Y.dat" ), "w" )

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
                write( handle_C2_X, string( distances[ index[ neigh+1] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ carbon, index[2] ] # dC-O1
            b = distance_matrix[ carbon, index[3] ] # dC-O2
            c = distance_matrix[ index[2], index[3] ] # dO1-O2
            # Angle through Al-Kashi
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_C2_X, string( alkash_angle, " ", step, "\n" ) )
            if check
                write( handle_C2_Y, string( distances[ index[ 2 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_angle_C2_Y, string( alkash_angle, " ", step, "\n" ) )
            end
        elseif nb_neighbors == 3
            check = true
            for neigh = 1:3
                write( handle_C3_X, string( distances[ index[ neigh+1] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
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
            alpha1 = geom.angleAlKash( a, b, c )
            write( handle_angle_C2, string( alpha1, " ", step, "\n" ) )
            write( handle_angle_C3_shortlong_X , string( alpha1, " ", step, "\n" ) )
            write( handle_angle_C3_23_X, string( alpha1, " ", step, "\n" ) )
            alpha2 = geom.angleAlKash( a, d, f )
            write( handle_angle_C2, string( alpha2, " ", step, "\n" ) )
            write( handle_angle_C3_shortlong_X , string( alpha2, " ", step, "\n" ) )
            write( handle_angle_C3_24_X, string( alpha2, " ", step, "\n" ) )
            alpha3 = geom.angleAlKash( b, d, e )
            write( handle_angle_C2, string( alpha3, " ", step, "\n" ) )
            write( handle_angle_C3_longlong_X , string( alpha3, " ", step, "\n" ) )
            write( handle_angle_C3_34_X, string( alpha3, " ", step, "\n" ) )
            #
            write( handle_C3_short_X, string( distances[ index[ 2 ] ], " " ) ) # neigh +1 because we ignore the 0
            if index[2] < nbC
                write( handle_C3_short_X, string( 0, " ", step, "\n" ) ) # neigh +1 because we ignore the 0
            else
                write( handle_C3_short_X, string( 1, " ", step, "\n" ) ) # neigh +1 because we ignore the 0
            end
            write( handle_C3_long_X, string( distances[ index[ 3 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
            write( handle_C3_long_X, string( distances[ index[ 4 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
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
            normal_v = normal_v./norm( normal_v)
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist = abs( dot( v_oc, normal_v ) )
            write( handle_base_C3_X, string( dist, " ", step, "\n") ) # Writting to disk
            if check
                write( handle_C3_Y, string( distances[ index[ 2 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_Y, string( distances[ index[ 3 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_Y, string( distances[ index[ 4 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_short_Y, string( distances[ index[ 2 ] ] , " " ) ) # neigh +1 because we ignore the 0
                if index[2] < nbC
                    write( handle_C3_short_Y, string( 0, " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                else
                    write( handle_C3_short_Y, string( 1, " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                end
                write( handle_C3_long_Y, string( distances[ index[ 3 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_C3_long_Y, string( distances[ index[ 4 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C3_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C3_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C3_Y, string( alkash_angle, " ", step, "\n" ) )
                write( handle_angle_C3_shortlong_Y , string( alpha1, " ", step, "\n" ) )
                write( handle_angle_C3_shortlong_Y , string( alpha2, " ", step, "\n" ) )
                write( handle_angle_C3_longlong_Y , string( alpha3, " ", step, "\n" ) )
                write( handle_angle_C3_23_X, string( alpha1, " ", step, "\n" ) )
                write( handle_angle_C3_24_X, string( alpha2, " ", step, "\n" ) )
                write( handle_angle_C3_34_X, string( alpha3, " ", step, "\n" ) )
                write( handle_base_C3_Y, string( dist, " ", step, "\n") )
            end
        elseif nb_neighbors >= 4

            check = true
            for neigh = 1:4
                write( handle_C4_X, string( distances[ index[ neigh+1] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end

            write( handle_C4_short_X, string( distances[2], " ", step, "\n" ) )
            write( handle_C4_long_X, string( distances[3], " ", step, "\n" ) )
            write( handle_C4_long_X, string( distances[4], " ", step, "\n" ) )

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

            angle23 = geom.angleAlKash( a, b, c )
            write( handle_angle_C4, string( angle23, " ", step, "\n" ) )
            write( handle_angle_C4_23_X, string( angle23, " ", step, "\n" ) )

            angle24 = geom.angleAlKash( a, d, f )
            write( handle_angle_C4, string( angle24, " ", step, "\n" ) )
            write( handle_angle_C4_24_X, string( angle24, " ", step, "\n" ) )

            angle34 = geom.angleAlKash( b, d, e )
            write( handle_angle_C4, string( angle34, " ", step, "\n" ) )
            write( handle_angle_C4_34_X, string( angle34, " ", step, "\n" ) )

            angle25 = geom.angleAlKash( a, g, j )
            write( handle_angle_C4, string( angle25, " ", step, "\n" ) )
            write( handle_angle_C4_34_X, string( angle25, " ", step, "\n" ) )

            angle35 = geom.angleAlKash( b, g, h )
            write( handle_angle_C4, string( angle35, " ", step, "\n" ) )
            write( handle_angle_C4_35_X, string( angle35, " ", step, "\n" ) )

            angle45 = geom.angleAlKash( d, g, i )
            write( handle_angle_C4, string( angle45, " ", step, "\n" ) )
            write( handle_angle_C4_45_X, string( angle45, " ", step, "\n" ) )

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
            normal_v = normal_v./norm( normal_v)
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist1 = abs( dot( v_oc, normal_v ) )
            write( handle_base_C4_X, string( dist1, " ", step, "\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 3, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            normal_v /= norm( normal_v)
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist2 = abs( dot( v_oc, normal_v ) )
            write( handle_base_C4_X, string( dist2, " ", step, "\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 2, : ] - neigh_positions[ 3, : ]
            v2 = neigh_positions[ 2, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            normal_v = normal_v./norm( normal_v)
            # Compute vector C-O1
            v_oc = neigh_positions[ 2, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist3 = abs( dot( v_oc, normal_v ) )
            write( handle_base_C4_X, string( dist3, " ", step, "\n") ) # Writting to disk
            #
            v1 = neigh_positions[ 1, : ] - neigh_positions[ 2, : ]
            v2 = neigh_positions[ 1, : ] - neigh_positions[ 4, : ]
            # Computing the normal to the plan
            normal_v = cross( v1, v2 )
            normal_v = normal_v./norm( normal_v)
            # Compute vector C-O1
            v_oc = neigh_positions[ 1, : ] - center_C
            # The distance to the base is the scalar product of C-O1 by the norm
            dist4 = abs( dot( v_oc, normal_v ) )
            write( handle_base_C4_X, string( dist4, " ", step, "\n") ) # Writting to disk
            if check
                write( handle_C4_Y, string( distances[ index[ 2 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 3 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 4 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_C4_Y, string( distances[ index[ 5 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( a, g, j )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( b, g, h )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                alkash_angle = geom.angleAlKash( d, g, i )
                write( handle_angle_C4_Y, string( alkash_angle, " ", step, "\n" ) )
                write( handle_angle_C4_23_Y, string( angle23, " ", step, "\n" ) )
                write( handle_angle_C4_24_Y, string( angle24, " ", step, "\n" ) )
                write( handle_angle_C4_25_Y, string( angle25, " ", step, "\n" ) )
                write( handle_angle_C4_34_Y, string( angle34, " ", step, "\n" ) )
                write( handle_angle_C4_35_Y, string( angle35, " ", step, "\n" ) )
                write( handle_angle_C4_45_Y, string( angle45, " ", step, "\n" ) )
                write( handle_base_C4_Y, string( dist1, " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist2, " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist3, " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_base_C4_Y, string( dist4, " ", step, "\n") ) # neigh +1 because we ignore the 0
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
                write( handle_O1_X, string( distances[ index[ neigh+1] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            if check
                write( handle_O1_Y, string( distances[ index[ neigh+1] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
            end
        elseif nb_neighbors == 2
            check = true
            for neigh = 1:2
                write( handle_O2_X, string( distances[ index[ neigh+1] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                if distances[ index[ neigh + 1 ]  ] > cut_off_low
                    check = false
                end
            end
            a = distance_matrix[ nbC+oxygen, index[2] ] # dO-C1
            b = distance_matrix[ nbC+oxygen, index[3] ] # dO-C2
            c = distance_matrix[ index[2], index[3] ] # dC1-C2
            alkash_angle = geom.angleAlKash( a, b, c )
            write( handle_angle_O2, string( alkash_angle, " ", step,  "\n" ) )
            if check
                write( handle_O2_Y, string( distances[ index[ 2 ] ], " ", step, "\n" ) ) # neigh +1 because we ignore the 0
                write( handle_O2_Y, string( distances[ index[ 3 ] ], " ", step, "\n") ) # neigh +1 because we ignore the 0
                write( handle_angle_O2_Y, string( alkash_angle, " ", step, "\n" ) )
            end
        else
            continue
        end
    end
    contact_matrix.writeStepMatrix( handle_matrix, distance_matrix )
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

close( handle_angle_C3_23_X )
close( handle_angle_C3_24_X )
close( handle_angle_C3_34_X )

close( handle_angle_C3_23_Y )
close( handle_angle_C3_24_Y )
close( handle_angle_C3_34_Y )

close( handle_angle_C4 )
close( handle_angle_C4_Y )

close( handle_angle_C4_23_X )
close( handle_angle_C4_24_X )
close( handle_angle_C4_25_X )
close( handle_angle_C4_34_X )
close( handle_angle_C4_35_X )
close( handle_angle_C4_45_X )

close( handle_angle_C4_23_Y )
close( handle_angle_C4_24_Y )
close( handle_angle_C4_25_Y )
close( handle_angle_C4_34_Y )
close( handle_angle_C4_35_Y )
close( handle_angle_C4_45_Y )

close( handle_angle_O2 )
close( handle_angle_O2_Y )

close( handle_angle_C3_shortlong_X )
close( handle_angle_C3_longlong_X )

close( handle_angle_C3_shortlong_Y )
close( handle_angle_C3_longlong_Y )

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



nb_box = 50
max_step = 20000
block_size = 100
min_nb = 5
max_angle = 179.5

data = readData( string( folder_out, "C2_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C2_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "C3_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C3_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "C4_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C4_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "C2_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C2_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "C3_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C3_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "C4_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "C4_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "O1_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "O1_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "O2_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "O2_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end


data = readData( string( folder_out, "O1_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "O1_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "O2_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "O2_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "distance_C3_short_X.dat" ), 3 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C3_short_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C3_long_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C3_long_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "distance_C3_short_Y.dat" ), 3 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C3_short_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C3_long_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C3_long_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C4_short_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C4_short_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C4_long_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C4_long_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C4_short_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C4_short_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "distance_C4_long_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "distance_C4_long_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "base_dist_C3_X.dat"  ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "base_C3_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "base_dist_C3_Y.dat"  ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "base_C3_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "base_dist_C4_X.dat"  ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "base_C4_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "base_dist_C4_Y.dat"  ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    file_hist = string( folder_out, "base_C4_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC2_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    nb_box = size(hist_avg)[1]
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC2_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "angleC2_Y.dat" ), 1 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC2_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "angleC3_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "angleC4_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleO2_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleO2_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
data = readData( string( folder_out, "angleO2_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleO2_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end


data = readData( string( folder_out, "angleC3_X_shortlong.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_X_shortlong_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_X_longlong.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_X_longlong_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_Y_shortlong.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_Y_shortlong_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_Y_longlong.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_Y_longlong_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_23_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_23_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_24_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_24_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_34_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_34_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_23_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_23_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_24_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_24_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC3_34_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC3_34_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_23_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_23_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_24_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_24_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_25_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_25_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_34_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_34_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_35_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_35_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_45_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_45_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

handle_angle_C4_23_Y = open( string( folder_out, "angleC4_23_Y.dat" ), "w" )
data = readData( string( folder_out, "angleC4_45_X.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_45_X_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_24_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_24_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_25_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_25_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_34_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_34_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_35_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_35_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end

data = readData( string( folder_out, "angleC4_45_Y.dat" ), 2 )
if data != false && size(data)[1] > min_nb*block_size
    hist_avg, hist_std, delta_hist, min_hist = makeHist( data, nb_box, max_step, block_size )
    sum_ = 0
    for i_box=1:nb_box
        if i_box*delta_hist+min_hist < max_angle
            hist_avg[i_box] = hist_avg[i_box]/( sin( (i_box*delta_hist+min_hist)*conversion.degr2rad ) )
            global sum_ += hist_avg[i_box]
        end
    end
    for i_box=1:nb_box
        hist_avg[i_box] = hist_avg[i_box]/sum_
    end
    file_hist = string( folder_out, "angleC4_45_Y_hist.dat" )
    writeHist( file_hist, hist_avg, hist_std, delta_hist, min_hist )
end
