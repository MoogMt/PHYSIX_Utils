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

frame_num=50

path_to_file = string( path_to_trimer_folder, "ELF_", frame_num, ".cube")
output_path_file = string( path_to_trimer_folder, "aELF_", frame_num, "_treated_mod.cube" )

atoms, cell, volume = cube_mod.readCube( path_to_file )

# Carbon index
carbons = [ 7, 10, 14 ]

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
move_vector_real = center_ - ( carbon_barycenter - volume.origin )

# Converts move vector as Int, to move a given number of voxels
# - Initialize vector for move
move_vector = zeros(Int, 3 )
# Loop over dimensions
for i=1:3
    # Compute move vector
    move_vector[i] = round(Int, move_vector_real[i]/ ( norm( volume.vox_vec[:,i] )/volume.nb_vox[i] ) )
end

volume = cube_mod.moveValues( volume, move_vector )

carbon_position_grid = zeros(Int,3,3)

for atom=1:3
    carbon_position_grid[:,atom] = cube_mod.getClosestIndex( atoms.positions[:,atom], volume , cell )
end

cut_off = 1.5

volume = cube_mod.carveCube( volume, carbon_position_grid, cut_off )

cube_mod.writeCube( output_path_file, atoms, cell, volume )
