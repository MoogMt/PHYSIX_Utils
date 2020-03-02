# Compute the number of occurences of trimers
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix
using geom
using utils

max_step=20000

Volumes=[ 9.325, 9.35, 9.375, 9.4, 9.5, 9.8, 10.0 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

cut_off = 1.75

min_angle=0
max_angle=181
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
        if nb_lines == 0
            continue
        end

        print("V: ", V, " T: ", T, "K\n")

        index_C = zeros(Int, nb_lines, 3 )
        index_O = zeros(Int, nb_lines, 3 )
        steps = zeros(Int, nb_lines )

        handle_in = open( target_file )
        for line=1:nb_lines
            keywords = split( readline( handle_in ) )
            steps[ line ] = parse(Int, keywords[1] )
            for carbon=1:3
                index_C[ line, carbon ] = parse(Int, keywords[ carbon + 1 ] )
            end
            for oxygen=1:3
                index_O[ line, oxygen ] = parse(Int, keywords[ oxygen + 4 ] )
            end
        end
        close( handle_in )

        cell = cell_mod.Cell_param( V, V, V )

        traj=filexyz.readFileAtomList( target_traj )

        # Angle through Al-Kashi
        target2_file = string( folder_target, "trimer-",cut_off,"-anglesCOC_mol.dat" )
        target3_file = string( folder_target, "trimer-",cut_off,"-anglesOCC_mol.dat" )
        handle_out_2 = open( target2_file, "w" )
        handle_out_3 = open( target3_file, "w" )
        for step_=1:nb_lines
            # Angle O-C-O
            distC1O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 1 ] )
            distC1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 3 ] )
            distO1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 3 ] )
            angle_C1 = geom.angleAlKash( distC1O1, distC1O3, distO1O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C1, "\n" ) )
            anglesOCO[ round( Int, ( angle_C1 - min_angle )/delta_angle ) + 1 ] += 1

            distC2O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 1 ] )
            distC2O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 2 ] )
            distO1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 2 ] )
            angle_C2 = geom.angleAlKash( distC2O1, distC2O2, distO1O2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C2, "\n" ) )
            anglesOCO[ round( Int, ( angle_C2 - min_angle )/delta_angle ) + 1 ] += 1

            distC3O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 2 ] )
            distC3O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 3 ] )
            distO2O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 2 ], index_O[ step_, 3 ] )
            angle_C3 = geom.angleAlKash( distC3O2, distC3O3, distO2O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C3, "\n" ) )
            anglesOCO[ round( Int, ( angle_C3 - min_angle )/delta_angle ) + 1 ] += 1

            # Angle C-O-C
            distC1C2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 2 ] )
            distC1C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 3 ] )
            distC2C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_C[ step_, 3 ] )

            angle_O1 = geom.angleAlKash( distC1O1, distC2O1, distC1C2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O1, "\n" ) )
            anglesCOC[ round( Int, (angle_O1 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O2 = geom.angleAlKash( distC2O2, distC3O2, distC2C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O2, "\n" ) )
            anglesCOC[ round( Int, (angle_O2 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O3 = geom.angleAlKash( distC1O3, distC3O3, distC1C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O3, "\n" ) )
            anglesCOC[ round( Int, (angle_O3 - min_angle)/delta_angle ) + 1 ] += 1
        end
        close( handle_out_2 )
        close( handle_out_3 )
    end
end

for ibox = 5:nb_box-1
    anglesCOC[ ibox ] = anglesCOC[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
    anglesOCO[ ibox ] = anglesOCO[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
end
anglesCOC /= sum( anglesCOC )
anglesOCO /= sum( anglesOCO )

file_out_COC = open( string( folder_out, "anglesCOC_moleculars_trim.dat" ), "w" )
file_out_OCO = open( string( folder_out, "anglesOCO_moleculars_trim.dat" ), "w" )
for ibox = 1:nb_box
    Base.write( file_out_COC, string( ibox*delta_angle + min_angle, " ", anglesCOC[ ibox ], "\n" ) )
    Base.write( file_out_OCO, string( ibox*delta_angle + min_angle, " ", anglesOCO[ ibox ], "\n" ) )
end
close( file_out_COC )
close( file_out_OCO )


Volumes=[ 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

cut_off = 1.75

min_angle=0
max_angle=181
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
        if nb_lines == 0
            continue
        end

        index_C = zeros(Int, nb_lines, 3 )
        index_O = zeros(Int, nb_lines, 3 )
        steps = zeros(Int, nb_lines )

        handle_in = open( target_file )
        for line=1:nb_lines
            keywords = split( readline( handle_in ) )
            steps[ line ] = parse(Int, keywords[1] )
            for carbon=1:3
                index_C[ line, carbon ] = parse(Int, keywords[ carbon + 1 ] )
            end
            for oxygen=1:3
                index_O[ line, oxygen ] = parse(Int, keywords[ oxygen + 4 ] )
            end
        end
        close( handle_in )

        cell = cell_mod.Cell_param( V, V, V )

        traj=filexyz.readFileAtomList( target_traj )

        # Angle through Al-Kashi
        target2_file = string( folder_target, "trimer-",cut_off,"-anglesCOC_mol.dat" )
        target3_file = string( folder_target, "trimer-",cut_off,"-anglesOCC_mol.dat" )
        handle_out_2 = open( target2_file, "w" )
        handle_out_3 = open( target3_file, "w" )
        for step_=1:nb_lines
            # Angle O-C-O
            distC1O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 1 ] )
            distC1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 3 ] )
            distO1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 3 ] )
            angle_C1 = geom.angleAlKash( distC1O1, distC1O3, distO1O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C1, "\n" ) )
            anglesOCO[ round( Int, ( angle_C1 - min_angle )/delta_angle ) + 1 ] += 1

            distC2O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 1 ] )
            distC2O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 2 ] )
            distO1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 2 ] )
            angle_C2 = geom.angleAlKash( distC2O1, distC2O2, distO1O2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C2, "\n" ) )
            anglesOCO[ round( Int, ( angle_C2 - min_angle )/delta_angle ) + 1 ] += 1

            distC3O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 2 ] )
            distC3O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 3 ] )
            distO2O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 2 ], index_O[ step_, 3 ] )
            angle_C3 = geom.angleAlKash( distC3O2, distC3O3, distO2O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C3, "\n" ) )
            anglesOCO[ round( Int, ( angle_C3 - min_angle )/delta_angle ) + 1 ] += 1

            # Angle C-O-C
            distC1C2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 2 ] )
            distC1C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 3 ] )
            distC2C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_C[ step_, 3 ] )

            angle_O1 = geom.angleAlKash( distC1O1, distC2O1, distC1C2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O1, "\n" ) )
            anglesCOC[ round( Int, (angle_O1 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O2 = geom.angleAlKash( distC2O2, distC3O2, distC2C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O2, "\n" ) )
            anglesCOC[ round( Int, (angle_O2 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O3 = geom.angleAlKash( distC1O3, distC3O3, distC1C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O3, "\n" ) )
            anglesCOC[ round( Int, (angle_O3 - min_angle)/delta_angle ) + 1 ] += 1
        end
        close( handle_out_2 )
        close( handle_out_3 )
    end
end

for ibox = 5:nb_box-1
    anglesCOC[ ibox ] = anglesCOC[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
    anglesOCO[ ibox ] = anglesOCO[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
end
anglesCOC /= sum( anglesCOC )
anglesOCO /= sum( anglesOCO )

file_out_COC = open( string( folder_out, "anglesCOC_polymeric_trim.dat" ), "w" )
file_out_OCO = open( string( folder_out, "anglesOCO_polymeric_trim.dat" ), "w" )
for ibox = 1:nb_box
    Base.write( file_out_COC, string( ibox*delta_angle + min_angle, " ", anglesCOC[ ibox ], "\n" ) )
    Base.write( file_out_OCO, string( ibox*delta_angle + min_angle, " ", anglesOCO[ ibox ], "\n" ) )
end
close( file_out_COC )
close( file_out_OCO )


Volumes=[ 9.325, 9.35, 9.375, 9.4, 9.5, 9.8, 10.0, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]

cut_off = 1.75

min_angle=0
max_angle=181
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
        if nb_lines == 0
            continue
        end

        print("V: ", V, " T: ", T, "K\n")

        index_C = zeros(Int, nb_lines, 3 )
        index_O = zeros(Int, nb_lines, 3 )
        steps = zeros(Int, nb_lines )

        handle_in = open( target_file )
        for line=1:nb_lines
            keywords = split( readline( handle_in ) )
            steps[ line ] = parse(Int, keywords[1] )
            for carbon=1:3
                index_C[ line, carbon ] = parse(Int, keywords[ carbon + 1 ] )
            end
            for oxygen=1:3
                index_O[ line, oxygen ] = parse(Int, keywords[ oxygen + 4 ] )
            end
        end
        close( handle_in )

        cell = cell_mod.Cell_param( V, V, V )

        traj=filexyz.readFileAtomList( target_traj )

        # Angle through Al-Kashi
        target2_file = string( folder_target, "trimer-",cut_off,"-anglesCOC_mol.dat" )
        target3_file = string( folder_target, "trimer-",cut_off,"-anglesOCC_mol.dat" )
        handle_out_2 = open( target2_file, "w" )
        handle_out_3 = open( target3_file, "w" )
        for step_=1:nb_lines
            # Angle O-C-O
            distC1O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 1 ] )
            distC1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_O[ step_, 3 ] )
            distO1O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 3 ] )
            angle_C1 = geom.angleAlKash( distC1O1, distC1O3, distO1O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C1, "\n" ) )
            anglesOCO[ round( Int, ( angle_C1 - min_angle )/delta_angle ) + 1 ] += 1

            distC2O1 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 1 ] )
            distC2O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_O[ step_, 2 ] )
            distO1O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 1 ], index_O[ step_, 2 ] )
            angle_C2 = geom.angleAlKash( distC2O1, distC2O2, distO1O2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C2, "\n" ) )
            anglesOCO[ round( Int, ( angle_C2 - min_angle )/delta_angle ) + 1 ] += 1

            distC3O2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 2 ] )
            distC3O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 3 ], index_O[ step_, 3 ] )
            distO2O3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_O[ step_, 2 ], index_O[ step_, 3 ] )
            angle_C3 = geom.angleAlKash( distC3O2, distC3O3, distO2O3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_C3, "\n" ) )
            anglesOCO[ round( Int, ( angle_C3 - min_angle )/delta_angle ) + 1 ] += 1

            # Angle C-O-C
            distC1C2 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 2 ] )
            distC1C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 1 ], index_C[ step_, 3 ] )
            distC2C3 = cell_mod.distance( traj[ steps[ step_ ] ], cell, index_C[ step_, 2 ], index_C[ step_, 3 ] )

            angle_O1 = geom.angleAlKash( distC1O1, distC2O1, distC1C2 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O1, "\n" ) )
            anglesCOC[ round( Int, (angle_O1 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O2 = geom.angleAlKash( distC2O2, distC3O2, distC2C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O2, "\n" ) )
            anglesCOC[ round( Int, (angle_O2 - min_angle)/delta_angle ) + 1 ] += 1

            angle_O3 = geom.angleAlKash( distC1O3, distC3O3, distC1C3 )
            Base.write( handle_out_2, string( steps[step_], " ", angle_O3, "\n" ) )
            anglesCOC[ round( Int, (angle_O3 - min_angle)/delta_angle ) + 1 ] += 1
        end
        close( handle_out_2 )
        close( handle_out_3 )
    end
end

for ibox = 5:nb_box-1
    anglesCOC[ ibox ] = anglesCOC[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
    anglesOCO[ ibox ] = anglesOCO[ ibox ]/( sind( (ibox-0.5)*delta_angle + min_angle ) )
end
anglesCOC /= sum( anglesCOC )
anglesOCO /= sum( anglesOCO )

file_out_COC = open( string( folder_out, "anglesCOC_all_trim.dat" ), "w" )
file_out_OCO = open( string( folder_out, "anglesOCO_all_trim.dat" ), "w" )
for ibox = 1:nb_box
    Base.write( file_out_COC, string( ibox*delta_angle + min_angle, " ", anglesCOC[ ibox ], "\n" ) )
    Base.write( file_out_OCO, string( ibox*delta_angle + min_angle, " ", anglesOCO[ ibox ], "\n" ) )
end
close( file_out_COC )
close( file_out_OCO )
