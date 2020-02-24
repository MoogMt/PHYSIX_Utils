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

folder_out = string( folder_target, "Data/Geometry" )
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

handle_matrix = open( string( folder_out_matrix, "distance_matrix.dat" ), "w" )

nb_step = size(traj)[1]
nb_atoms=sum(species_nb)

cell = cell_mod.Cell_param(V,V,V)


for step=1:nb_step
    distance_matrix = contact_matrix.buildMatrix( traj[step], cell)
    for carbon=1:nbC
        distances = distance_matrix[carbon,:]
        nb_neighbors = sum(ones(size(distances)[1])[ distances .< cut_off ] ) - 1
        index_nearest = sortperm(distances)
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
            write( handle_angle_C2, string( alkash_angle, " " ) )
            if check
                write( handle_C2_Y, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                write( handle_angle_C2_Y, string( alkash_angle, " " ) )
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
            write( handle_angle_C2, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( a, d, f )
            write( handle_angle_C2, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( b, d, e )
            write( handle_angle_C2, string( alkash_angle, " " ) )
            if check
                write( handle_C3_Y, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C3_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C3_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C3_Y, string( alkash_angle, " " ) )
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
            write( handle_angle_C4, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( a, d, f )
            write( handle_angle_C4, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( b, d, e )
            write( handle_angle_C4, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( a, g, j )
            write( handle_angle_C4, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( b, g, h )
            write( handle_angle_C4, string( alkash_angle, " " ) )
            alkash_angle = geom.angleAlKash( d, g, i )
            write( handle_angle_C4, string( alkash_angle, " " ) )
            if check
                write( handle_C4_Y, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
                alkash_angle = geom.angleAlKash( a, b, c )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( a, d, f )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( b, d, e )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( a, g, j )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( b, g, h )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
                alkash_angle = geom.angleAlKash( d, g, i )
                write( handle_angle_C4_Y, string( alkash_angle, " " ) )
            end
        end
    end
    for oxygen=1:nbO
        distances = distance_matrix[ nbC+oxygen,:]
        nb_neighbors = sum( ones( size( distances )[1] )[ distances .< cut_off ] ) - 1
        index_nearest = sortperm( distances )
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
            if check
                write( handle_O2_Y, string( distances[ index[ neigh+1] ] ,"\n") ) # neigh +1 because we ignore the 0
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

close( handle_matrix )

# end

function getCO3( file_out, nb_steps, nb_atoms , nbC, nbO,  nb_steps, )
    fileC=open( file_out ,"w")
    # Loop over steps
    for step=1:nb_steps
        print("Progres: ",step/nb_steps*100,"%\n")
        # Loop over carbons
        for carbon=1:nbC
            # Keep distances of oxygen to carbon
            distances = zeros(nbO)
            for oxygen=1:nbO
                distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
            end
            # get the index of the oxygen sorted by distances
            index=sortperm(distances)

            # write distances to file
            write(fileC, string(step*unit," ") )
            for i=1:max_neigh
                write(fileC, string(distances[index[i]]," ") )
            end

            # Computes and write the angles between the N nearest oxygen
            for i=1:max_neigh-1
                for j=i+1:max_neigh
                    a=cell_mod.distance(traj[step],cell,carbon,Int(index[i]+nbC))
                    b=cell_mod.distance(traj[step],cell,carbon,Int(index[j]+nbC))
                    c=cell_mod.distance(traj[step],cell,Int(index[i]+nbC),Int(index[j]+nbC))
                    write(fileC,string(acosd((a*a+b*b-c*c)/(2*a*b))," "))
                end
            end

            write(fileC,string("\n"))
        end
    end
    close(fileC)
end
