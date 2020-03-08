# Using code from LibAtomicSim
using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod

# Using code from standard Julia Libraries
using Statistics

# Computes the coordinances for all files

function countCoordinances( matrix::Array{T1,2}, max_number_neighbor ) where {T1 <: Real }
    nb_atoms=size(matrix)[1]
    count_coord = zeros(Int, max_number_neighbor+1)
    for atom=1:nb_atoms-1
        count_coord[ sum( matrix[atom,:] )+1  ] += 1
    end
    return count_coord
end
function countCoordinances( file_path::T1, V::T2, cut_off::T3, max_number_neighbor::T4  ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Int }
    traj = filexyz.readFileAtomList( file_path )
    if traj == false
        return 0, 0, 0, 0
    end
    cell = cell_mod.Cell_param(V,V,V)
    nb_step = atom_mod.getNbStep( traj )
    species  = atom_mod.getSpecies( traj[1] )
    nb_species = atom_mod.getNbElementSpecies( traj[1], species )
    start_species = atom_mod.getStartSpecies( traj[1], species )
    n_spec = size(species)[1]
    count_coordinances=zeros(Int, nb_step, max_number_neighbor+1, n_spec )
    for step=1:nb_step
        distance_matrix=contact_matrix.buildMatrix( traj[step], cell,cut_off )
        for i_spec=1:n_spec
            start_ = start_species[i_spec]
            end_ = start_species[i_spec]+nb_species[i_spec]-1
            count_coordinances[ step, : , i_spec ] = countCoordinances( distance_matrix[start_:end_,:], max_number_neighbor )
        end
    end
    return count_coordinances, species, nb_species, start_species
end
function writeRawCounts( prefix_out::T1, counts::Array{T2,3}, species::Vector{T3}, cut_off::T4  ) where { T1 <: AbstractString, T2 <: Int, T3 <: AbstractString, T4 <: Real }
    n_spec = size(species)[1]
    nb_step = size(counts)[1]
    n_max  = size(counts)[2]
    for i_spec = 1:n_spec
        handle_out = open( string( prefix_out, "_", species[i_spec], "_", cut_off, ".dat"), "w" )
        for step=1:nb_step
            Base.write( handle_out, string( step, " " ) )
            for n_=1:n_max
                Base.write( handle_out, string( counts[step,n_,i_spec], " " ) )
            end
            Base.write( handle_out, string( "\n" ) )
        end
        close( handle_out )
    end
    return true
end
function writeRawFracs( prefix_out::T1, counts::Array{T2,3}, species::Vector{T3}, nb_species::Vector{T4}, cut_off::T5  ) where { T1 <: AbstractString, T2 <: Int, T3 <: AbstractString, T4 <: Int, T5 <: Real }
    n_spec = size(species)[1]
    nb_step = size(counts)[1]
    n_max  = size(counts)[2]
    for i_spec = 1:n_spec
        handle_out = open( string( prefix_out, "_", species[i_spec], "_", cut_off, ".dat"), "w" )
        for step=1:nb_step
            Base.write( handle_out, string( step, " " ) )
            for n_=1:n_max
                Base.write( handle_out, string( round( counts[step,n_,i_spec]/nb_species[i_spec]*100, digits=3 ), " " ) )
            end
            Base.write( handle_out, string( "\n" ) )
        end
        close( handle_out )
    end
    return true
end
function countRawCoordinances( folder_target::T1, prefix_out_count::T2, prefix_out_frac::T3, max_nb_neighbor::T4, cut_off::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: Int, T5 <: Real }

    traj_file = string(folder_target, "TRAJEC_fdb_wrapped.xyz")
    if ! isfile( traj_file )
        return false
    end

    counts, species, nb_species, start_species  = countCoordinances( traj_file , V, cut_off, max_nb_neighbor  )
    if counts == false
        return false
    end

    # Output is Data folder
    folder_out_loc=string( folder_target, "Data/Coordinance/")
    if ! isdir( folder_out_loc )
        Base.Filesystem.mkdir( folder_out_loc )
    end

    # Computing Raw Counts
    prefix_out_raw=string(folder_out_loc,prefix_out_count)
    writeRawCounts( prefix_out_raw, counts, species, cut_off )

    # Computing Fractions
    prefix_out_frac_raw = string(folder_out_loc,prefix_out_frac )
    writeRawFracs( prefix_out_frac_raw, counts, species, nb_species, cut_off )

    return true
end
function readCount( file_target::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_target )
        return false
    end
    handle_in = open( file_target)
    # Count the number of step and max neighbor
    nb_step = 0
    nb_max = 0
    while ! eof( handle_in )
        test=readline(handle_in)
        if nb_step == 0
            nb_max = size( split( test ) )[1]-1
        end
        nb_step+=1
    end
    # Read
    seekstart(handle_in) # back at the begining
    counts_raw = zeros(nb_step,nb_max)
    for step=1:nb_step
        keyword = split( readline( handle_in) )
        for i=1:nb_max
            counts_raw[step,i] = parse(Int, keyword[i+1] )
        end
    end
    close(handle_in)
    return counts_raw
end
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


cut_off = 1.75
max_nb_neighbor = 6
species=["C","O"]
species_nb=[32,64]
prefix_count_raw = string("counts_coord_raw")
prefix_frac_raw = string("frac_coord_raw")
folder_base=string("/media/mathieu/Elements/CO2/")


Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperatures=[3000,2500,2000]

# Computing Raw costs
restart_count = false
for T in Temperatures
    for V in Volumes
        folder_target=string(folder_base, V, "/",T,"K/")
        if ! isfile( string( folder_target, "Data/Coordinance/", prefix_count_raw, "_", species[1], "_",cut_off,".dat" ) ) || restart_count
            countRawCoordinances( folder_target, prefix_count_raw, prefix_frac_raw, max_nb_neighbor, cut_off )
        end
    end
end


# Computing statistics
n_spec = size(species)[1]
handle_out_temp=Array{IO}( undef, n_spec, 2 )
folder_out_general = string( folder_base, "Data/Coordinance/" )
folder_in_pressure = string( folder_base, "Data/Pressure/" )
if ! isdir( folder_out_general )
    Base.Filesystem.mkdir( folder_out_general )
end
for T in Temperatures
    for i_spec=1:n_spec
        handle_out_temp[i_spec,1] = open( string( folder_out_general, T, "K_coord_", species[i_spec], "stat_glob.dat" ), "w" )
        handle_out_temp[i_spec,2] = open( string( folder_out_general, T, "K_coord_", species[i_spec], "stat_div.dat" ), "w" )
    end
    pressure, delta_press, volumes = readPressure( string(folder_in_pressure,T,"_P_global.dat"))
    for V in Volumes
        folder_target=string(folder_base, V, "/",T,"K/")
        for i_spec=1:n_spec
            folder_out_loc = string(folder_target, "Data/Coordinance/")
            target_file = string( folder_out_loc, prefix_count_raw, "_", species[i_spec], "_",cut_off,".dat" )
            if ! isfile( target_file )
                continue
            end
            press = pressure[ volumes.== V ][1]
            dp = delta_press[ volumes.== V ][1]
            counts_ = readCount( target_file )
            max_nb = size(counts_)[2]
            nb_step = size(counts_)[1]
            coord_global=copy(counts_)
            coord_div=copy(counts_)
            for step=1:nb_step
                for i=1:max_nb
                    coord_global[step,i] = counts_[step,i]*(i-1)/species_nb[i_spec]
                    coord_div[step,i] = counts_[step,i]/species_nb[i_spec]
                end
            end
            handle_out_loc = open( string(folder_out_loc,"stat_global_coord_",species[i_spec],"_",cut_off,".dat"), "w" )
            Base.write( handle_out_temp[i_spec],  string( V, " ", Statistics.mean(coord_global)," ",Statistics.std(coord_global),"\n" ) )
            Base.write( handle_out_loc, string( Statistics.mean(coord_global)," ",Statistics.std(coord_global) ) )
            close( handle_out_loc )
            handle_out_loc_div = open( string(folder_out_loc,"stat_div_coord_",species[i_spec],"_",cut_off,".dat"), "w" )
            Base.write( handle_out_temp[i_spec,2], string( V, " " ))
            for i=1:max_nb
                Base.write( handle_out_temp[i_spec,2],  string( Statistics.mean(coord_div[:,i])," ",Statistics.std(coord_div[:,i])," " ) )
                Base.write( handle_out_loc_div, string( Statistics.mean(coord_div[:,i])," ",Statistics.std(coord_div[:,i])," " ) )
            end
            Base.write( handle_out_temp[i_spec,2], string("\n" ) )
            Base.write( handle_out_loc_div, string( "\n" ) )
            close( handle_out_loc_div )
        end
    end
    for i_spec=1:size(species)[1]
        close( handle_out_temp[i_spec,1] )
        close( handle_out_temp[i_spec,2] )
    end
end
