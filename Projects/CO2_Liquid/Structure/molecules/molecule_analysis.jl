# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
using graph
using exp_data
using geom

function writeMolecule( handle_out::T1 , positions::Array{T2}, species::Vector{T3}, step::T4, mol_index::T5 ) where { T1 <: IO, T2 <: Real,  T3 <: AbstractString, T4 <: Int, T5 <: Int }
    size_molecule=size(positions)[1]
    Base.write( handle_out, string( size_molecule, "\n" ) )
    Base.write( handle_out, string( "STEP: ",step," ", mol_index, "\n" ) )
    for atom=1:size_molecule
        Base.write( handle_out, string(species[atom], " ") )
        for i=1:3
            Base.write( handle_out, string( round(positions[atom,i],digits=3), " " ) )
        end
        Base.write( handle_out, string("\n") )
    end
    return true
end
function writeSizeHistogramTime( file_out::T1, hist_time::Array{T2,2} )  where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w" )
    nb_step = size(hist_time)[1]
    nb_atoms  = size(hist_time)[2]
    for step = 1:nb_step
        Base.write( handle_out, string( step,  " " ) )
        for size_atom=1:nb_atoms
            Base.write( handle_out, string( hist_time[ step, size_atom ], " " ) )
        end
        Base.write( handle_out, string("\n") )
    end
    close( handle_out )
end
function writeSizeHistogramGlobal( file_out::T1, hist_time::Array{T2,2} )  where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w" )
    nb_atoms  = size(hist_time)[2]
    for size_atom=1:nb_atoms
        Base.write( handle_out, string( size_atom, " ", sum(hist_time[ :, size_atom ]), "\n" ) )
    end
    close( handle_out )
end

# Thermodynamical values
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperature=[1750,2000,2500,3000,3500]
cut_off=1.75

V=8.82
T=3000


folder_base = string( "/media/mathieu/Elements/CO2/" )
folder_target = string( folder_base, V, "/", T, "K/" )

file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )

traj = filexyz.readFileAtomList( file_traj )
cell = cell_mod.Cell_param(V,V,V)
cut_off=1.75


nb_atoms = size(traj[1].names)[1]
nb_step = size(traj)[1]

folder_target_mol = string( folder_target, "Data/Molecules/" )

if ! isdir( folder_target_mol )
    Base.Filesystem.mkdir( folder_target_mol )
end

handle_mol_inf = Vector{ IO }( undef, nb_atoms )
handle_mol_fin = Vector{ IO }( undef, nb_atoms )
for atom=1:nb_atoms
  handle_mol_inf[atom]= open( string( folder_target_mol, "mol_inf_",atom,".xyz"), "w" )
  handle_mol_fin[atom]= open( string( folder_target_mol, "mol_fin_",atom,".xyz"), "w" )
end
# handle_mol_inf_link = open( string( folder_target_mol, "mol_inf_link.xyz" ) , "w" )

nb_step  = size( traj )[1]
nb_atoms = size( traj[1].names )[1]

hist_fin=zeros( Int, nb_step, nb_atoms )
hist_inf=zeros( Int, nb_step, nb_atoms )
hist_gen=zeros( Int, nb_step, nb_atoms )

for step=1:nb_step

    print("Progress: ",round(step/nb_step*100,digits=3),"%\n")

    # Computing the molecules through graph exploration
    positions_local=copy(traj[step].positions)
    matrix = contact_matrix.buildMatrix( traj[step], cell, cut_off )
    molecules=graph.getGroupsFromMatrix(matrix)
    nb_molecules=size(molecules)[1]
    matrices=graph.extractAllMatrixForTrees( matrix, molecules )

    # Loop over mmolecules
    for molecule=1:nb_molecules
        size_molecule = size( molecules[molecule] )[1]
        # Molecule is an atom, discard
        if size_molecule <= 1
            continue
        end
        # Check if molecule is finite
        visited=zeros(Int,size(molecules[molecule]))
        adjacent_molecule=getAllAdjacentVertex(matrices[molecule])
        #isinf, isok = cell_mod.isInfiniteChain( visited, matrices[molecule], adjacent_molecule, positions_local, cell, 1, molecules[molecule], cut_off )
        isinf = cell_mod.checkInfiniteChain( matrices[molecule], positions_local, cell, molecules[molecule], cut_off )
        # If something went wrong, moving on...
        # if ! isok
        #     print( "Issue reconstructing molecule ", molecule, " at step: ", step, "\n" )
        #     continue
        # end
        # Whether infinite or finite, the size of the molecule is counted, and it is printed in a specific file
        if isinf
            writeMolecule( handle_mol_inf[size_molecule], positions_local[molecules[molecule],:], traj[step].names[ molecules[molecule] ], step, molecule ) # write molecule
            # If the molecule is infinite, we print the atoms that should be bonded but aren't.
            #list = cell_mod.findUnlinked( matrices[molecule], positions_local[molecules[molecule],:], cell, molecules[molecule], cut_off )
            # We put the positions of those atoms in a specific file, for visualization purposes (with VMD)
            #writeMolecule( handle_mol_inf_link, positions_local[unique(list),:], traj[step].names[ molecules[unique(list)] ], step, molecule ) # write links
            # Counting sizes
            hist_inf[ step, size_molecule ] += 1
        else
            writeMolecule( handle_mol_fin[size_molecule], positions_local[molecules[molecule],:], traj[step].names[ molecules[molecule] ], step, molecule ) # write molecule
            # Counting sizes
            hist_fin[ step, size_molecule ] += 1
        end
        hist_gen[ step, size_molecule ] += 1
        # NB: here we only count the molecule by size, however that isn't fair, as a large molecule may take the
        # totality of the atoms. Therefore, a corrective factor will be needed when computing fair statistics.

    end
end

for atom=1:nb_atoms
    close( handle_mol_inf[atom] )
    close( handle_mol_fin[atom] )
    if hist_fin[ atom ] == 0
        Base.Filesystem.rm( string( folder_target_mol, "mol_fin_",atom,".xyz") )
    end
    if hist_inf[ atom ] == 0
        Base.Filesystem.rm( string( folder_target_mol, "mol_inf_",atom,".xyz") )
    end
end

writeSizeHistogramTime( string( folder_target_mol, "hist_inf.dat"), hist_inf )
writeSizeHistogramTime( string( folder_target_mol, "hist_fin.dat"), hist_fin )
writeSizeHistogramTime( string( folder_target_mol, "hist_gen.dat"), hist_gen )

writeSizeHistogramGlobal( string( folder_target_mol, "hist_inf.dat"), hist_inf )
writeSizeHistogramGlobal( string( folder_target_mol, "hist_fin.dat"), hist_fin )
writeSizeHistogramGlobal( string( folder_target_mol, "hist_gen.dat"), hist_gen )
