computer="win"

if computer == "win"
    push!(LOAD_PATH,"C:\\Users\\moogm\\OneDrive\\Documents\\GitHub\\LibAtomicSim\\Julia\\")
end

using utils
using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using cube_mod

function lineBetweenAtoms( atoms::T1, cell_params::T2, volume::T3, nb_points::T4, index_atom1::T5, index_atom2::T6 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: cube_mod.Volume, T4 <: Int, T5 <: Int, T6 <: Int }
    # Arguments
    # - atoms: AtomList (from atom_mod) containing atomic positions, index and names
    # - cell_params: Cell_param (from cell_mod) containing cell parameters
    # - volume: Volume (from cube_mod) containing all ELF/Density information
    # - nb_points: number of points to take in the line between the atoms
    # - index_atom1, index_atom2: indexes of the target atoms
    # Output
    # - distFrom1: distance from first point
    # - data: data on the line

    # Initialize vector for 1 to 2 displacement
    delta_vector_from_1to2 = zeros(Real, 3)

    # Loop over dimensions
    for i=1:3
        # Compute vector from 1 to 2 s
        delta_vector_from_1to2[i] = atoms.positions[i, index_atom2] - atoms.positions[i, index_atom1]

        # Apply PBC
        if delta_vector_from_1to2[i] > cell_params.length[i]*0.5
            delta_vector_from_1to2[i] = delta_vector_from_1to2[i] - cell_params.length[i]
        end
        if delta_vector_from_1to2[i] < - cell_params.length[i]*0.5
            delta_vector_from_1to2[i] = delta_vector_from_1to2[i] + cell_params.length[i]
        end
    end

    # Compute displacement vector to go from 1 to 2
    delta_vector_from_1to2 = delta_vector_from_1to2 / nb_points

    # Initialize vector for output data
    data = zeros(Real, nb_points)
    # Loop over points
    for i=1:nb_points
        # Get the position, wrapping point if need be
        positions_line = cell_mod.wrap( atoms.positions[:, index_atom1] + i*delta_vector_from_1to2 - volume.origin, cell_params )

        # Get the indexs of the closests grid point for target position
        indexs = getClosestIndex( positions_line, volume )

        # Get the data at the desired point
        data[i] = volume.matrix[ indexs[1], indexs[2], indexs[3] ]
    end

    return data
end

folder_base=""

if computer == "mac"
    folder_base=string("/Volumes/DATA2/PHYSIX/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/ELF/8.82/Trajectory_2/")
elseif computer == "win"
    folder_base=string("F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\8.82\\Trajectory_2\\")
elseif computer == "ubuntu"
    folder_base=string("/media/moogmt/DATA2/PHYSIX/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/ELF/8.82/Trajectory_2/")
end

# Number of structures to loop over
max_structure=1000

# Number of points to take along a C-O line
nb_points_lines=50

# Number of atoms
nbC=32
nbO=64

# Max number of neighbors to consider around carbon
max_neigh=5

# Parameters for bonding
cut_off_bonded     = 1.6
cut_off_not_bonded = 1.9

# Initialize data tensors
data_distances = zeros(Real, max_neigh, nbC, max_structure )
data_bonded    = zeros(Real, max_neigh, nbC, max_structure )
data_dens      = zeros(Real, nb_points_lines, max_neigh, nbC, max_structure )
data_elf       = zeros(Real, nb_points_lines, max_neigh, nbC, max_structure )

# Writting data to file
file_out_input_dens = string( folder_base, "input_dens.dat" )
file_out_input_elf  = string( folder_base, "input_elf.dat"  )

# Opening output files
handle_input_dens_out  = open( file_out_input_dens, "w" )
handle_input_elf_out   = open( file_out_input_elf,  "w" )

# Computing distance
for step_ = 1:max_structure
    # Write progress
    print("Progress: ", step_/max_structure*100, "%\n" )

    # Create path for the input path for density and elf files
    dens_path = string( folder_base, step_, "_density.cube" )
    elf_path  = string( folder_base, step_, "_elf.cube"     )

    # Read ELF and DENSITY data
    atoms, cell_params, volume_dens = cube_mod.readCubeVolume( dens_path )
    atoms, cell_params, volume_elf  = cube_mod.readCubeVolume( elf_path  )

    # Initialize local data repository
    nearest_neighbors = zeros(Real, max_neigh, nbC )
    dist_nearest      = zeros(Real, max_neigh, nbC )

    # Loop over carbons
    for carbon=1:nbC
        # Local storing data
        distance_oxygen = zeros(Real, nbO )

        # Loop over oxygen atoms
        for oxygen=1:nbO
            distance_oxygen[oxygen] = cell_mod.distanceOrtho( atoms.positions[:,carbon], atoms.positions[:,nbC+oxygen], cell_params.length )
        end

        # Compute neighbors of local carbon
        nearest_neighbors[:, carbon] = sortperm( distance_oxygen )[1:max_neigh] .+ nbC

        # Loop over neighbors
        for neigh=1:max_neigh
            # Get ELF and DENSITY values
            data_dens[ :, neigh, carbon, step_ ] = lineBetweenAtoms( atoms, cell_params, volume_dens, nb_points_lines, carbon, nearest_neighbors[ neigh, carbon ] )
            data_elf[ :, neigh, carbon, step_ ]  = lineBetweenAtoms( atoms, cell_params, volume_elf,  nb_points_lines, carbon, nearest_neighbors[ neigh, carbon ] )

            # Computing distances between carbon and its neighbor
            dist =  cell_mod.distanceOrtho( atoms.positions[ :, carbon ], atoms.positions[ :, nearest_neighbors[ neigh, carbon ] ], cell_params.length )

            # Storing distance
            data_distances[ neigh, carbon, step_ ] = dist

            # Determining if bonded or not or unknown
            if dist < cut_off_bonded
                data_bonded[ neigh, carbon, step_ ] = 1
            elseif dist > cut_off_not_bonded
                data_bonded[ neigh, carbon, step_ ] = -1
            end

            # Begining of line
            line_start = string( step_, " ", carbon, " ", neigh, " ", nearest_neighbors[neigh, carbon], " ", dist, " " )

            # Write basic carbon information
            write( handle_input_dens_out, line_start )
            write( handle_input_elf_out,  line_start )

            # - Loop over points
            for point=1:nb_points_lines
                # Writting data to files
                write( handle_input_dens_out, string( data_dens[ point, neigh, carbon, step_ ], " " ) )
                write( handle_input_elf_out,  string( data_elf[ point, neigh, carbon, step_ ],  " " ) )
            end

            # Writting output data to file
            write( handle_input_dens_out, string( data_bonded[ neigh, carbon, step_ ], "\n" ) )
            write( handle_input_elf_out,  string( data_bonded[ neigh, carbon, step_ ], "\n" ) )
        end
    end
end

# Closing files
close( handle_input_dens_out )
close( handle_input_elf_out  )

# Histogram preparation
nb_box = 100

# Determines histogram parameters for density
min_value_dens = minimum(data_dens)
max_value_dens = maximum(data_dens)
delta_box_dens = ( max_value_dens - min_value_dens )/nb_box

# Determines histogram parameters for ELF
min_value_elf = minimum( data_elf )
max_value_elf = maximum( data_elf )
delta_box_elf = ( max_value_elf - min_value_elf )/nb_box

# Initialize histograms tables for density
hist_all_dens       = zeros(Real, nb_box + 1, nb_points_lines )
hist_neighbors_dens = zeros(Real, nb_box + 1, nb_points_lines, max_neigh )
hist_bonded_dens    = zeros(Real, nb_box + 1, nb_points_lines, 3 )

# Initialize histograms tables for ELF
hist_all_elf       = zeros(Real, nb_box + 1, nb_points_lines )
hist_neighbors_elf = zeros(Real, nb_box + 1, nb_points_lines, max_neigh )
hist_bonded_elf    = zeros(Real, nb_box + 1, nb_points_lines, 3 )

# Making Histogram for both density and ELF
for point=1:nb_points_lines
    for neigh=1:max_neigh
        for carbon=1:nbC
            for structure=1:max_structure
                #
                hist_all_dens[       round(Int, ( data_dens[ point, neigh, carbon, structure ] - min_value_dens )/delta_box_dens ) + 1, point ]        += 1
                hist_neighbors_dens[ round(Int, ( data_dens[ point, neigh, carbon, structure ] - min_value_dens )/delta_box_dens ) + 1, point, neigh ] += 1
                #
                hist_all_elf[       round(Int, ( data_elf[ point, neigh, carbon, structure ] - min_value_elf )/delta_box_elf ) + 1, point ]        += 1
                hist_neighbors_elf[ round(Int, ( data_elf[ point, neigh, carbon, structure ] - min_value_elf )/delta_box_elf ) + 1, point, neigh ] += 1

                if data_distances[ neigh, carbon, structure ] < cut_off_bonded
                    hist_bonded_dens[ round(Int, ( data_dens[ point, neigh, carbon, structure ] - min_value_dens )/delta_box_dens ) + 1,  point, 1 ] += 1
                    hist_bonded_elf[ round(Int,  (  data_elf[ point, neigh, carbon, structure ] - min_value_elf )/delta_box_elf )   + 1,  point, 1 ] += 1
                elseif data_distances[ neigh, carbon, structure ] > cut_off_not_bonded
                    hist_bonded_dens[ round(Int, ( data_dens[ point, neigh, carbon, structure ] - min_value_dens )/delta_box_dens ) + 1,  point, 2 ] += 1
                    hist_bonded_elf[ round(Int,  (  data_elf[ point, neigh, carbon, structure ] - min_value_elf )/delta_box_elf )   + 1,  point, 2 ] += 1
                else
                    hist_bonded_dens[ round(Int, ( data_dens[ point, neigh, carbon, structure ] - min_value_dens )/delta_box_dens ) + 1,  point, 3 ] += 1
                    hist_bonded_elf[ round(Int,  (  data_elf[ point, neigh, carbon, structure ] - min_value_elf )/delta_box_elf )   + 1,  point, 3 ] += 1
                end
            end
        end
        hist_neighbors_dens[ :, point, neigh ] = hist_neighbors_dens[ :, point, neigh ]/sum( hist_neighbors_dens[ :, point, neigh ] )
        hist_neighbors_elf[  :, point, neigh ] = hist_neighbors_elf[  :, point, neigh ]/sum( hist_neighbors_elf[  :, point, neigh ] )
    end
    hist_all_dens[ :, point ] = hist_all_dens[ :, point ]/sum( hist_all_dens[ :, point ] )
    hist_all_elf[ :, point ]  = hist_all_elf[  :, point ]/sum( hist_all_elf[  :, point ] )
    # Can be replaced by a loop
    hist_bonded_dens[ :, point, 1 ] = hist_bonded_dens[ :, point, 1 ]/sum( hist_bonded_dens[ :, point, 1 ] )
    hist_bonded_elf[  :, point, 1 ] = hist_bonded_elf[  :, point, 1 ]/sum( hist_bonded_elf[  :, point, 1 ] )
    hist_bonded_dens[ :, point, 2 ] = hist_bonded_dens[ :, point, 2 ]/sum( hist_bonded_dens[ :, point, 2 ] )
    hist_bonded_elf[  :, point, 2 ] = hist_bonded_elf[  :, point, 2 ]/sum( hist_bonded_elf[  :, point, 2 ] )
    hist_bonded_dens[ :, point, 3 ] = hist_bonded_dens[ :, point, 3 ]/sum( hist_bonded_dens[ :, point, 3 ] )
    hist_bonded_elf[  :, point, 3 ] = hist_bonded_elf[  :, point, 3 ]/sum( hist_bonded_elf[  :, point, 3 ] )
end

# Determining histogram files names
file_hist_all_dens = string( folder_base, "all_hist_density.dat" )
file_hist_all_elf  = string( folder_base, "all_hist_elf.dat" )
file_hist_neighbors_dens = string( folder_base, "neigbors_hist_density.dat" )
file_hist_neighbors_elf  = string( folder_base, "neigbors_hist_elf.dat" )
file_hist_bonded_dens = string( folder_base, "bonded_hist_dens.dat")
file_hist_bonded_elf = string( folder_base, "bonded_hist_elf.dat")

# Writting data out
# Opening output files
handle_out_dens  = open( file_hist_all_dens, "w" )
handle_out_elf   = open( file_hist_all_elf,  "w" )
handle_out2_dens = open( file_hist_neighbors_dens, "w" )
handle_out2_elf  = open( file_hist_neighbors_elf,  "w" )
handle_out3_dens = open( file_hist_bonded_dens, "w" )
handle_out3_elf  = open( file_hist_bonded_elf,  "w" )
# Loop over points
for point=1:nb_points_lines
    for box=1:nb_box
        write( handle_out_dens,  string( point, " ", round( box*delta_box_dens + min_value_dens, digits=3 ), " ", hist_all_dens[ box, point ], "\n" ) )
        write( handle_out2_dens, string( point, " ", round( box*delta_box_dens + min_value_dens, digits=3 ), " " ) )
        write( handle_out3_dens, string( point, " ", round( box*delta_box_dens + min_value_dens, digits=3 ), " " ) )
        write( handle_out_elf,   string( point, " ", round( box*delta_box_elf  + min_value_elf,  digits=3 ), " ", hist_all_elf[ box, point ], "\n" ) )
        write( handle_out2_elf,  string( point, " ", round( box*delta_box_elf  + min_value_elf,  digits=3 ), " " ) )
        write( handle_out3_elf,  string( point, " ", round( box*delta_box_elf  + min_value_elf,  digits=3 ), " " ) )
        for neigh=1:max_neigh
            write( handle_out2_dens, string( hist_neighbors_dens[ box, point, neigh ], " " ) )
            write( handle_out2_elf,  string( hist_neighbors_elf[ box, point, neigh ],  " " ) )
        end
        for i=1:3
            write( handle_out3_dens, string( hist_bonded_dens[ box, point, i ], " " ) )
            write( handle_out3_elf,  string( hist_bonded_elf[ box, point, i ],  " " ) )
        end
        write( handle_out2_dens, string("\n") )
        write( handle_out2_elf, string("\n") )
        write( handle_out3_dens, string("\n") )
        write( handle_out3_elf, string("\n") )
    end
    write( handle_out_dens, "\n")
    write( handle_out2_dens, string("\n") )
    write( handle_out3_dens, string("\n") )
    write( handle_out_elf, "\n")
    write( handle_out2_elf, string("\n") )
    write( handle_out3_elf, string("\n") )
end
# Closing files
close( handle_out_dens  )
close( handle_out2_dens )
close( handle_out3_dens )
close( handle_out_elf   )
close( handle_out2_elf  )
close( handle_out3_elf  )
