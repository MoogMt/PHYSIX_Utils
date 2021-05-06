path_to_libatomsim = "C:\\Users\\moogm\\OneDrive\\Documents\\GitHub\\LibAtomicSim\\Julia\\"
push!( LOAD_PATH, path_to_libatomsim )

using atom_mod
using cell_mod
using geom
using cube_mod
using periodicTable
using conversion

using LinearAlgebra

path_to_trimer_folder = "F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\Rings\\"

frames = [9, 50, 86, 91, 96]

for frame_num in frames

    print( "Handling frame: ", frame_num, "\n" )
    path_to_file = string( path_to_trimer_folder, "ELF_", frame_num, ".cube")
    output_path_file = string( path_to_trimer_folder, "aELF_", frame_num, "_treated_mod.cube" )
    output_path_file_2 = string( path_to_trimer_folder, "bELF_", frame_num, "_treated_mod.cube" )
    output_path_file_3 = string( path_to_trimer_folder, "cELF_", frame_num, "_treated_mod.cube" )

    # Reading volume
    atoms, cell, volume = cube_mod.readCube( path_to_file )

    # Rewrite the same volume
    cube_mod.writeCube( output_path_file_3, atoms, cell, volume )

    # Index of the carbon in the trimer
    carbons = [ 7, 10, 14 ]

    # Index of the oxygen in the trimer
    oxygens = [ 72, 60, 59, 51, 52, 46 ]

    # Compute position of the barycenter of carbon atoms
    # NB not perfect, should account for PBC first, but that will do for now
    # Initialize vector for carbon barycenter
    carbon_barycenter = zeros(Real, 3 )
    # - Loop over dimensions
    for i=1:3
        # Loop over carbon atoms
        for carbon in carbons
            carbon_barycenter[i] += atoms.positions[ i, carbon ]
        end
        # Average
        carbon_barycenter[i] /= 3
    end

    # Compute center of the cell
    # NB Not perfect, should improve to acccount for non-orthorombic cells
    # - Initialize vector for center of cell
    center_ = zeros(Real, 3 )
    # Loop over dimensions
    for i=1:3
        center_[i] = cell.length[i]/2
    end

    # Compute the move vector
    move_vector_real = center_ - ( carbon_barycenter + volume.origin )

    # Moving atoms in order to recenter
    # - Loop over atoms
    for atom=1:size(atoms.positions)[2]
        # Loop over dimensions
        for i=1:3
            # Moving atoms in coordinate i
            atoms.positions[i,atom] = atoms.positions[i,atom] + move_vector_real[i]
        end
    end

    # Converts move vector as Int, to move a given number of voxels
    # - Initialize vector for move
    move_vector = zeros(Int, 3 )
    # Loop over dimensions
    for i=1:3
        # Compute move vector
        move_vector[i] = round(Int, move_vector_real[i]/ ( norm( volume.vox_vec[:,i] )/volume.nb_vox[i] ) )
    end

    # Recentring values of the grid
    volume = cube_mod.moveValues( volume, move_vector )

    cube_mod.writeCube( output_path_file_2, atoms, cell, volume )

    atoms_positions_grid = zeros(Int,3,9)
    # Computing grid positions of target atoms
    # - Adding the position in the grid of carbon atom
    for i=1:3
        atoms_positions_grid[:,i]  = cube_mod.getClosestIndex( atoms.positions[:, carbons[i] ], volume )
    end
    # - Adding the piosition in the grid of oxygen atoms
    for i=1:6
        atoms_positions_grid[:,i+3] = cube_mod.getClosestIndex( atoms.positions[:, oxygens[i] ], volume )
    end

    # Determining cut-off for carving
    cut_off = 1.0
    fac = 1.1
    soft_fac = 0.5

    # Carving data within a given cut-off of atoms
    volume = cube_mod.carveCube( volume, atoms_positions_grid, cut_off, fac, soft_fac )

    # Writting new volume to file
    cube_mod.writeCube( output_path_file, atoms, cell, volume )
end
