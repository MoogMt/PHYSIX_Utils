module cube_mod

include("atoms.jl");
include("cell.jl")

# Reads a cube file and returns all or parts of its informations
function readCube{T1<:AbstractString}( file_name::T1)
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    # Number of atoms
    nb_atoms = parse(Int,split(lines[3])[1]);

    # Oigin position of the density
    #-----------------------------------------------------
    center=Vector{Real}(3)
    for i=1:3
        center[i] = parse(Float64, split( lines[3] )[i+1] );
    end
    #-----------------------------------------------------

    # Number of voxels in each direction
    #-----------------------------------------------------
    nb_vox=Vector{Real}(3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[3+i] )[1] )
    end
    #-----------------------------------------------------

    # Reads Cell Matrix
    #-----------------------------------------------------
    cell_matrix=cell_mod.Cell_matrix();
    for i=1:3
        for j=1:3
            cell_matrix.matrix[i,j] = parse(Float64, split( lines[4+i] )[1+j] )*0.52917721067
        end
    end
    #-----------------------------------------------------

    # Reads atoms
    atom_list=AtomList();
    

    # Reads Density

    return nb_atoms, nb_vox, cell_matrix
end

end
