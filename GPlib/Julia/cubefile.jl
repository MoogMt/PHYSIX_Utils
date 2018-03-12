module cube_mod

include("atoms.jl");
include("cell.jl");

#-----------------------------------------------------------------------
mutable struct Volume
    position::Array{Real}
    value::Vector{Real}
    function Volume()
        new(Array{Real}(0,3),Vector{Real}())
    end
    function Volume{T2 <: Real}( nb_vox::Vector{T2} )
        if size(nb_vox)[1] == 3
            nb_tot=nb_vox[1]*nb_vox[2]*nb_vox[3]
            new(Array{Real}(nb_tot,3),Vector{Real}(nb_tot))
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
end
#-----------------------------------------------------------------------

# Reads a cube file and returns all or parts of its informations
function readCube{T1<:AbstractString}( file_name::T1)
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    #------------------
    # Number of atoms
    #-----------------------------------------
    nb_atoms = parse(Int,split(lines[3])[1]);
    #-----------------------------------------

    #--------------------------------
    # Oigin position of the density
    #-----------------------------------------------------
    center=Vector{Real}(3)
    for i=1:3
        center[i] = parse(Float64, split( lines[3] )[i+1] );
    end
    #-----------------------------------------------------

    #-------------------------------------
    # Number of voxels in each direction
    #-----------------------------------------------------
    nb_vox=Vector{Real}(3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[3+i] )[1] )
    end
    #-----------------------------------------------------

    #--------------------
    # Reads Cell Matrix
    #-----------------------------------------------------
    cell_matrix=cell_mod.Cell_matrix();
    for i=1:3
        for j=1:3
            cell_matrix.matrix[i,j] = parse(Float64, split( lines[4+i] )[1+j] )*0.52917721067
        end
    end
    #-----------------------------------------------------

    #-------------
    # Reads atoms
    #----------------------------------------------------
    atom_list=AtomList(nb_atoms);
    for i=1:nb_atoms
        atom_list.names[i] = split( lines[6+i] )[1]
        for j=1:3
            atom_list.position[i,j] = parse(Float64, split( lines[6+i])[2+j] )
        end
    end
    #----------------------------------------------------

    #----------------
    # Reads Density
    #----------------------------------------------------
    nb_tot=nb_vox[1]*nb_vox[2]*nb_vox[3]
    volume = Volume(nb_tot)
    for i=1:nb_tot/6
        for j=1:6
            
        end
    end
    #----------------------------------------------------

    return atom_list, cell_matrix, volume
end

end
