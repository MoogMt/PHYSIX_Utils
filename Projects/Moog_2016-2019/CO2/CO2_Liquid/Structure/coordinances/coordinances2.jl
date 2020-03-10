# Using code from LibAtomicSim
using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod

# Using code from standard Julia Libraries
using Statistics

function readPressure( file_pressure::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_pressure )
        return false, false, false
    end
    handle_press_in = open( file_pressure )
    nb_p = 0
    while !eof( handle_press_in )
        readline( handle_press_in )
        nb_p += 1
    end
    seekstart( handle_press_in )
    pressure = zeros(Real, nb_p )
    delta_p  = zeros(Real, nb_p )
    volumes  = zeros(Real, nb_p )
    for i=1:nb_p
        keyword=split(readline( handle_press_in ))
        pressure[i] = parse( Float64, keyword[3] )
        delta_p[i]  = parse( Float64, keyword[4] )
        volumes[i]  = parse( Float64, keyword[1] )
    end
    close( handle_press_in )
    return pressure, delta_p, volumes
end
function computeCoordinancesC( bond_matrix::Array{T1,2}  ) where { T1 <: Real }
    nbC=32
    coord = zeros(Int, nbC)
    for atom=1:nbC
        coord[atom] = Int( sum(bond_matrix[atom,:]) )
    end
    ones_dum = ones(Int, nbC )
    coord2 = sum( ones_dum[ coord .<= 2 ] )
    coord3 = sum( ones_dum[ coord .== 3 ] )
    coord4 = sum( ones_dum[ coord .>= 4 ] )
    return coord2, coord3, coord4
end
function computeCoordinancesO( bond_matrix::Array{T1,2} ) where { T1 <: Real }
    nbC=32
    nbO=64
    coord = zeros(Int, nbO)
    for atom=1:nbO
        coord[atom] = Int( sum( bond_matrix[nbC+atom,:]) )
    end
    ones_dum = ones(Int, nbO )
    coord1 = sum( ones_dum[ coord .== 1 ] )
    coord2 = sum( ones_dum[ coord .== 2 ] )
    return coord1, coord2
end


cut_off = 1.75
species=["C","O"]
species_nb=[32,64]

nb_step = 20000

fracC=100/32
fracO=100/64

Temperatures=[3000,2500,2000]

folder_base = string("/media/mathieu/Elements/CO2/")

for T in Temperatures
    file_pressure = string( folder_base, "Data/Pressure/",T,"_P_global.dat")
    P,dP,Volumes = readPressure( file_pressure )
    folder_out = string( folder_base, "Data/Coordinance/")
    if ! isdir( folder_out )
        Base.Filesystem.mkdir( folder_out )
    end
    file_out_C2 = open( string( folder_out, "C2_", T, ".dat" ), "w" )
    file_out_C3 = open( string( folder_out, "C3_", T, ".dat" ), "w" )
    file_out_C4 = open( string( folder_out, "C4_", T, ".dat" ), "w" )
    file_out_O1 = open( string( folder_out, "O1_", T, ".dat" ), "w" )
    file_out_O2 = open( string( folder_out, "O2_", T, ".dat" ), "w" )
    for i_vol=1:size(Volumes)[1]
        print("Volume: ",i_vol," ",Volumes[i_vol],"\n")
        cell = cell_mod.Cell_param( Volumes[i_vol], Volumes[i_vol], Volumes[i_vol] )
        folder_target = string( folder_base, Volumes[i_vol],"/",T,"K/")
        file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )
        if ! isfile( file_traj )
            continue
        end
        traj = filexyz.readFileAtomList( file_traj )
        C_coord = zeros( nb_step, 3 )
        O_coord = zeros( nb_step, 2 )
        for step=1:nb_step
            bm = contact_matrix.buildMatrix( traj[step], cell, cut_off )
            C_coord2, C_coord3, C_coord4 = computeCoordinancesC( bm )
            C_coord[step,1] = C_coord2*fracC
            C_coord[step,2] = C_coord3*fracC
            C_coord[step,3] = C_coord4*fracC
            O_coord1, O_coord2 = computeCoordinancesO( bm )
            O_coord[step,1] = O_coord1*fracO
            O_coord[step,2] = O_coord2*fracO
        end
        #
        C_coord2_avg = Statistics.mean( C_coord[:,1] )
        C_coord3_avg = Statistics.mean( C_coord[:,2] )
        C_coord4_avg = Statistics.mean( C_coord[:,3] )
        O_coord1_avg = Statistics.mean( O_coord[:,1] )
        O_coord2_avg = Statistics.mean( O_coord[:,2] )
        #
        C_coord2_std = Statistics.std( C_coord[:,1] )
        C_coord3_std = Statistics.std( C_coord[:,2] )
        C_coord4_std = Statistics.std( C_coord[:,3] )
        O_coord1_std = Statistics.std( O_coord[:,1] )
        O_coord2_std = Statistics.std( O_coord[:,2] )

        Base.write( file_out_C2, string( P[i_vol], " ", dP[i_vol], " ", C_coord2_avg, " ", C_coord2_std, "\n" ) )
        Base.write( file_out_C3, string( P[i_vol], " ", dP[i_vol], " ", C_coord3_avg, " ", C_coord3_std, "\n" ) )
        Base.write( file_out_C4, string( P[i_vol], " ", dP[i_vol], " ", C_coord4_avg, " ", C_coord4_std, "\n" ) )
        Base.write( file_out_O1, string( P[i_vol], " ", dP[i_vol], " ", O_coord1_avg, " ", O_coord1_std, "\n" ) )
        Base.write( file_out_O2, string( P[i_vol], " ", dP[i_vol], " ", O_coord2_avg, " ", O_coord2_std, "\n" ) )
    end
    close( file_out_C2 )
    close( file_out_C3 )
    close( file_out_C4 )
    close( file_out_O1 )
    close( file_out_O2 )
end
