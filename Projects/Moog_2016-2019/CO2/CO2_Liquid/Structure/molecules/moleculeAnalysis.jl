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

using LinearAlgebra

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
function getMaxDistance( positions::Array{T1,2} ) where { T1 <: Real }
    max_dist=0
    nb_atoms=size(positions)[1]
    for atom1=1:nb_atoms-1
        for atom2=atom1+1:nb_atoms
            distance=LinearAlgebra.norm( positions[atom1,:]-positions[atom2,:] )
            if distance > max_dist
                max_dist = distance
            end
        end
    end
    return max_dist
end
function getMoleculeLength( atoms::T1 ) where { T1 <: atom_mod.AtomList }
    return getMaxDistance( atoms.positions )
end
function getMoleculeLength( traj::Vector{T1} ) where { T1 <: atom_mod.AtomList }
    nb_step=size(traj)[1]
    lengths=Vector{Real}(undef,0)
    for step=1:nb_step
        push!( lengths, getMaxDistance( traj[step].positions) )
    end
    return lengths
end

# Thermodynamical values
Volumes=[ 10.0, 9.8, 9.5, 9.4, 9.375, 9.35, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6 ]
Temperatures=[ 1750, 2000, 2500, 3000 ]
cut_off=1.75


folder_base = string( "/media/mathieu/Elements/CO2/" )
folder_data_gen = string( folder_base, "/Data/Molecules/")
folder_in_pressure = string( folder_base, "Data/Pressure/" )
if  !isdir( folder_data_gen)
    Base.Filesystem.mkdir(folder_data_gen)
end
handle_inf = open( string(folder_data_gen,"molecules_inf_map.dat"),"w")
for T in Temperatures
    print("T: ",T,"\n")
    pressure, delta_press, volumes = readPressure( string(folder_in_pressure,T,"_P_global.dat"))
    for V in Volumes
        folder_target = string( folder_base, V, "/", T, "K/" )
        file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )
        traj = filexyz.readFileAtomList( file_traj )
        if traj == false
            continue
        end
        cell = cell_mod.Cell_param(V,V,V)
        nb_atoms = size(traj[1].names)[1]
        nb_step = size(traj)[1]
        folder_target_mol = string( folder_target, "Data/Molecules/" )
        if ! isdir( folder_target_mol )
            Base.Filesystem.mkdir( folder_target_mol )
        end
        #
        handle_mol_inf = Vector{ IO }( undef, nb_atoms )
        handle_mol_fin = Vector{ IO }( undef, nb_atoms )
        for atom=1:nb_atoms
            handle_mol_inf[atom]= open( string( folder_target_mol, "mol_inf_",atom,".xyz"), "w" )
            handle_mol_fin[atom]= open( string( folder_target_mol, "mol_fin_",atom,".xyz"), "w" )
        end
        #
        nb_step  = size( traj )[1]
        nb_atoms = size( traj[1].names )[1]
        # histograms of sizes
        hist_fin=zeros( Int, nb_step, nb_atoms )
        hist_inf=zeros( Int, nb_step, nb_atoms )
        hist_gen=zeros( Int, nb_step, nb_atoms )
        count_inf = 0
        # Analyzing molecules
        for step=1:nb_step
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
                # Whether infinite or finite, the size of the molecule is counted, and it is printed in a specific file
                if isinf
                    writeMolecule( handle_mol_inf[size_molecule], positions_local[molecules[molecule],:], traj[step].names[ molecules[molecule] ], step, molecule ) # write molecule
                    hist_inf[ step, size_molecule ] += 1
                else
                    writeMolecule( handle_mol_fin[size_molecule], positions_local[molecules[molecule],:], traj[step].names[ molecules[molecule] ], step, molecule ) # write molecule
                    hist_fin[ step, size_molecule ] += 1
                end
                hist_gen[ step, size_molecule ] += 1
                # NB: here we only count the molecule by size, however that isn't fair, as a large molecule may take the
                # totality of the atoms. Therefore, a corrective factor will be needed when computing fair statistics.
            end
            if sum( hist_inf[ step, : ] ) > 1
                count_inf += 1
            end
        end
        # Removing useless files
        for atom=1:nb_atoms
            close( handle_mol_inf[atom] )
            close( handle_mol_fin[atom] )
            if sum(hist_fin[ :,atom ]) == 0
                Base.Filesystem.rm( string( folder_target_mol, "mol_fin_",atom,".xyz") )
            end
            if sum(hist_inf[ :, atom ]) == 0
                Base.Filesystem.rm( string( folder_target_mol, "mol_inf_",atom,".xyz") )
            end
        end
        # Writting The dsitrbution of size at each step
        writeSizeHistogramTime( string( folder_target_mol, "hist_inf.dat"), hist_inf )
        writeSizeHistogramTime( string( folder_target_mol, "hist_fin.dat"), hist_fin )
        writeSizeHistogramTime( string( folder_target_mol, "hist_gen.dat"), hist_gen )
        # Writting the global size over all steps
        writeSizeHistogramGlobal( string( folder_target_mol, "hist_inf.dat"), hist_inf )
        writeSizeHistogramGlobal( string( folder_target_mol, "hist_fin.dat"), hist_fin )
        writeSizeHistogramGlobal( string( folder_target_mol, "hist_gen.dat"), hist_gen )
        # Computing Pressure
        press = pressure[ volumes.== V ][1]
        if sum(hist_inf[:,:]) > 0
            Base.write( handle_inf, string( press, " ", V, " ", T, " ", round( 1, digits=3 ), "\n" ) )
        else
            Base.write( handle_inf, string( press, " ", V, " ", T, " ", round( 0, digits=3 ), "\n" ) )
        end
    end
    write(handle_inf,string("\n"))
end
close(handle_inf)

# Computing length and comparing
nb_atoms=96
nb_box=50
for V in Volumes
    for T in Temperatures
        folder_target = string( folder_base, V, "/", T, "K/" )
        if ! isfile( string(folder_target, "TRAJEC_fdb.xyz") )
            continue
        end
        print("Processing: ",V," ",T,"K\n")
        folder_out = string( folder_target, "Data/Molecules/")
        lengths_total=Vector{Real}(undef,0)
        handle_out_avg = open( string( folder_out, "lengths_avg_mol_fin.dat"), "w")
        for size_ = 1:nb_atoms
            target_file = string( folder_out, "mol_fin_",size_,".xyz")
            if ! isfile( target_file )
                if isfile( string(folder_target_mol, "lengths_global_mol_fin_",size_,".dat") )
                    Base.Filesystem.rm( string(folder_target_mol, "lengths_global_mol_fin_",size_,".dat") )
                end
                if isfile( string(folder_target_mol, "hist_mol_fin_all.dat") )
                    Base.Filesystem.rm( string(folder_target_mol, "hist_mol_fin_all.dat") )
                end
                continue
            end
            molecules = filexyz.readFileAtomList( target_file )
            if molecules == false
                continue
            end
            lengths=getMoleculeLength(molecules)
            if typeof(molecules) == AtomList
                push!( lengths_total, lengths )
                Base.write( handle_out_avg, string( size_," ",lengths," ",0,"\n" ) )
                file_out_lengths = open( string(folder_target_mol, "lengths_mol_fin_",size_,".dat"), "w")
                Base.write( file_out_lengths, string(lengths, "\n"))
                close( file_out_lengths )
            else
                # Add the sizes to the general size vector
                for i=1:size(lengths)[1]
                    push!(lengths_total,lengths[i])
                end
                # Prints the average and standard deviation of lengths for a given size
                Base.write( handle_out_avg, string( size_, " ", Statistics.mean(lengths), " ", Statistics.std(lengths), "\n" ) )
                max_=maximum(lengths)
                min_=minimum(lengths)
                delta_=(max_-min_)/nb_box
                hist_nb=zeros(Real,nb_box+1)
                # Write all lengths for the size in a specific file and compute length histogram for size
                file_out_lengths = open( string(folder_target_mol, "lengths_mol_fin_",size_,".dat"), "w")
                for i=1:size(lengths)[1]
                    Base.write( file_out_lengths, string(lengths[i], "\n"))
                    hist_nb[ Int( trunc( (lengths[i]-min_)/delta_ )+1 ) ] += 1
                end
                close( file_out_lengths )
                # Prints the histogram of lengths for the target size
                file_out_hist = open( string(folder_target_mol, "hist_mol_fin_",size_,".dat"), "w" )
                for ibox=1:nb_box
                    Base.write( file_out_hist, string( (ibox*delta_)+min_," ", hist_nb[ibox],"\n" ) )
                end
            end
        end
        close(handle_out_avg)
        max_=maximum(lengths_total)
        min_=minimum(lengths_total)
        delta_=(max_-min_)/nb_box
        hist_nb=zeros(Real,nb_box+1)
        file_out_lengths = open( string(folder_target_mol, "lengths_global_mol_fin.dat"), "w")
        for i=1:size(lengths_total)[1]
            Base.write( file_out_lengths, string(lengths_total[i], "\n"))
            hist_nb[ Int( trunc( (lengths_total[i]-min_)/delta_ )+1 ) ] += 1
        end
        close(file_out_lengths)
        file_out_hist = open( string(folder_target_mol, "hist_mol_fin_all.dat"), "w" )
        for ibox=1:nb_box
            Base.write( file_out_hist, string( (ibox*delta_)+min_," ", hist_nb[ibox],"\n" ) )
        end
        close( file_out_hist )
    end
end
