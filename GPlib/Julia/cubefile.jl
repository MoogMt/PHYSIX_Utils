module cube_mod

include("atoms.jl");
include("cell.jl")

# Reads a cube file and returns all or parts of its informations
function readCube{T1<:AbstractString}( file_name::T1)
    # Open file
    file=open(file_name);
    lines=readlines(file);
    close(file);

    nb_atoms = parse(Int,split(lines[3])[1]);

    center=Vector{Real}(3)
    for i=1,3
        center[i] = parse(Real, split( lines[3] )[i+1] );
    end

    nb_vox=Vector{Real}(3)
    for i=1,3
        nb_vox[i] = parse(Real, split( lines[3+i] )[1] )
    end

    cell_matrix=Array{Real}(3,3)
    for i=1,3
        for j=1,3
            cell_matrix[i,j] = parse(Real, split( lines[4+i] )[1+j] )*0.52917721067
        end
    end
    cell_mod.Cell_matrix=cell_matrix

    

    return nb_atoms, nb_vox, cell_matrix
end

end
