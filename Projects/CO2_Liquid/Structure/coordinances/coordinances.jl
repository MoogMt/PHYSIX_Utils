GPfolder=string("/home/mathieu/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)

using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod

# Computes the coordinances for all files

function countCoordinances( matrix::Array{T1,2}, max_number_neighbor ) where {T1 <: Real }
    nb_atoms=size(matrix)[1]
    count_coord = zeros(Int, nb_atoms)
    for atom=1:nb_atoms-1
        count_coord[ sum( matrix[atom,:] ) + 1 ] += 1
    end
    return n_bonds
end

function countCoordinances( file_path::T1, V::T2, cut_off::T3, max_number_neighbor::T4  ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Int }
    traj = filexyz.readFileAtomList( file_path )
    if traj == false
        return 0, 0, 0, 0
    end
    cell = cell_mod.Cell_param(V,V,V)
    nb_step = atom_mod.getNbStep( traj )
    species  = atom_mod.getSpecies( traj[1] )
    nb_species = atom_mod.getNbSpecies( traj[1] )
    start_species = atom_mod.getStartSpecies( traj[1], species )
    count_coordinances=zeros( nb_step, max_number_neighbor+1, nb_species )
    for step=1:nb_step
        distance_matrix=contact_matrix.buildMatrix( traj[step], cell,cut_off )
        for i_spec=1:nb_species
            count_coordinances[ step, : , i_spec ] = countCoordinances( distance_matrix[ start_species[i_spec]:start_species[i_spec]+nb_species[i_spec], : ], max_number_neighbor )
        end
    end
    return count_coordinances, species, nb_species, start_species
end

cut_off = 1.75
max_nb_neighbor = 4

V=8.82
T=3000

folder_base=string("/media/mathieu/Elements/CO2/")

folder_target=string(folder_base, V, "/",T,"K/")

traj_file = string(folder_target, "TRAJEC_fdb_wrapped.xyz")

countCoordinances( traj_file , V, cut_off, max_nb_neighbor  )
