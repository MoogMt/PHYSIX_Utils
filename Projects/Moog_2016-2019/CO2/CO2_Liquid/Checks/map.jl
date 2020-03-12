using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using exp_data
using LsqFit
using Statistics
using exp_data


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

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base="/media/moogmt/Elements/CO2/"
folder_base="/media/mathieu/Elements/CO2/"

Volumes      = [ 10.0, 9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.1,9.05, 9.0, 8.82, 8.6 ]
Temperatures = [ 2000, 2500, 3000 ]

folder_out_map = string( folder_base, "Data/Map/" )
if ! isdir( folder_out_map )
    Base.Filesystem.mkdir( folder_out_map )
end

handle_map = open( string( folder_out_map,"map.dat" ), "w" )
for T in Temperatures
    file_pressure = string( folder_base, "Data/Pressure/",T,"_P_global.dat")
    P,dP,Volumes = readPressure( file_pressure )
    for i_vol=1:size(Volumes)[1]
        Base.write( handle_map, string( P[i_vol], " ", T,  "\n" ) )
    end
end
close( handle_map )
