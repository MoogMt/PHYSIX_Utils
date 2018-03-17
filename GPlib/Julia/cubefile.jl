include("cell.jl")
include("geom.jl")

module cube_mod

import atom_mod.AtomList
import cell_mod.Cell_matrix
importall atom_mod
importall cell_mod
importall geom

#-----------------------------------------------------------------------
mutable struct Volume
    matrix::Array{Real}
    function Volume()
        new( Array{Real}(0,4) )
    end
    function Volume{T1 <: Real}( nb_vox_iso::T1 )
        if nb_vox_iso > 0
            new( Array{Real}(nb_vox_iso*nb_vox_iso*nb_vox_iso,4))
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
    function Volume{T1 <: Real}( nb_vox::Vector{T1} )
        if size(nb_vox)[1] == 3
            nb_tot = nb_vox[1]*nb_vox[2]*nb_vox[3];
            new( Array{Real}(nb_tot,4) );
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
    # Origin position of the density
    #-----------------------------------------------------
    center=Vector{Real}(3)
    for i=1:3
        center[i]=parse(Float64,split(lines[3])[1+i])
    end
    #-----------------------------------------------------

    #-------------------------------------
    # Number of voxels in each direction
    #-----------------------------------------------------
    nb_vox=Vector{Int}(3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[3+i] )[1] )
    end
    #-----------------------------------------------------

    #--------------------
    # Reads Cell Matrix
    #-----------------------------------------------------
    cell_matrix=Cell_matrix();
    for i=1:3
        for j=1:3
            cell_matrix.matrix[i,j] = parse(Float64, split( lines[3+i] )[1+j] )*0.52918
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
            atom_list.positions[i,j] = parse(Float64, split( lines[6+i])[2+j] )
        end
    end
    #----------------------------------------------------

    #----------------
    # Reads Density
    #----------------------------------------------------
    nb_tot=nb_vox[1]*nb_vox[2]*nb_vox[3]
    nb_col=6
    volume = cube_mod.Volume(nb_vox)
    offset=(Int)(6+nb_atoms+1)
    x=0.; y=0.; z=0.;
    for i=0:nb_tot/nb_col-1
        for j=1:nb_col
            index=(Int)(i*nb_col+j)
            volume.matrix[index,1] = x
            volume.matrix[index,2] = y
            volume.matrix[index,3] = z
            volume.matrix[index,4] = parse(Float64, split( lines[(Int)(offset+i)])[j] )
            z=z+1;
            if z == nb_vox[3]
                z=0;
                y=y+1;
            end
            if y == nb_vox[2]
                y=0;
                z=0;
                x=x+1;
            end
        end
    end
    #----------------------------------------------------

    #------------------------------
    # Scaling positions for volume
    #----------------------------------------------------
    dV=[0.,0.,0.]
    for i=1:3
        for j=1:3
            dV[i] = dV[i] + cell_matrix.matrix[i,j];
        end
        volume.matrix[:,i] = volume.matrix[:,i]*dV[i]+center[i];
    end
    #----------------------------------------------------

    #-------------------------------------------
    # Scaling cell vectors by number of voxels
    #--------------------------------------------
    for i=1:3
        cell_matrix.matrix[i,i]=cell_matrix.matrix[i,i]*nb_vox[i]
    end
    #---------------------------------------------

    return atom_list, cell_matrix, volume
end

function traceVolume{ T1 <: Real, T2 <: Real, T3 <: Volume }( position1::Vector{T1}, position2::Vector{T2}, volume::T3 )
    if size(position1)[1] != 3 || size(position2)[1] != 3
        return false
    end
    indexs1= getClosest( position1, volume);
    indexs2= getClosest( position2, volume);
    dindex=index1-index2
    curseur=indexs1
    list=Array{Real}(0,3)
    while norm(curseur-indexs2) > 0
        vcat(index2,curseur)
        move(curseur,indexs2)
    end
end
print("Cube Module Loaded!\n")

end
