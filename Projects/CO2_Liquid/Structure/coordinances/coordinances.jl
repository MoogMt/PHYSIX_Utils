GPfolder=string("/home/mathieu/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)

using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod

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
function writeRawCounts( prefix_out::T1, counts::Array{T2,3}, species::Vector{T3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: AbstractString }
    n_spec = size(species)[1]
    nb_step = size(counts)[1]
    n_max  = size(counts)[2]
    for i_spec = 1:n_spec
        handle_out = open( string( prefix_out, "_", species[i_spec], ".dat"), "w" )
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

cut_off = 1.75
max_nb_neighbor = 5

V=8.82
T=3000

folder_base=string("/media/mathieu/Elements/CO2/")

folder_target=string(folder_base, V, "/",T,"K/")

traj_file = string(folder_target, "TRAJEC_fdb_wrapped.xyz")

counts, species, nb_species, start_species  = countCoordinances( traj_file , V, cut_off, max_nb_neighbor  )

folder_out_loc=string( folder_target, "Data/Coordinance/")

if ! isdir( folder_out_loc )
    Base.Filesystem.mkdir( folder_out_loc )
end

prefix_out_raw=string(folder_out_loc,"counts_coord_raw")
writeRawCounts( prefix_out_raw, counts, species )
