module cube_mod

include("atoms.jl");
include("cell.jl");

#-----------------------------------------------------------------------
mutable struct Volume
    matrix::Array{Real}
    function Volume()
        new( Array{Real}(0,4) )
    end
    function Volume{T2 <: Real}( nb_vox::Vector{T2} )
        if size(nb_vox)[1] == 3
            new( Array{Real}(nb_vox[1]*nb_vox[2]*nb_vox[3],4) )
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
dfvhsdefa    end
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
    nb_col=6
    volume = Volume(nb_tot)
    offset=6+nb_atoms
    x=0; y=0; z=0;
    for i=0:nb_tot/nb_col-1
        for j=1:nb_col
            volume.matrix[1,i*nb_col+j] = x
            volume.matrix[2,i*nb_col+j] = y
            volume.matrix[3,i*nb_col+j] = z
            volume.matrix[4,i*nb_col+j] = parse(Float64, split( lines[offset+i])[j] )
            z++
            if z == nb_atoms
                z=0;
                y=y+1;
            end
            if y == nb_atoms
                y=0;
                x=x+1;
            end
        end
    end
    # Scaling positions
    dV=[0,0,0]
    for i=1:3
        for j=1:3
            dV[i] += cell_matrix[i,j];
        end
        volume.matrix[i,:] = volume.matrix[i,:]*dV[i];
    end
    #----------------------------------------------------

    return atom_list, cell_matrix, volume
end

print("Cube Module Loaded!\n")

end
