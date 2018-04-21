include("cell.jl")
include("geom.jl")

module cube_mod

include("conversion.jl");

import atom_mod.AtomList
importall atom_mod
importall cell_mod
importall geom

#-----------------------------------------------------------------------
# Volume 4, nb_vox*nb_vox*nb_vox
#-----------------------
# 1-3 for positions
# 4 for the values
#-----------------------
mutable struct Volume
    matrix::Array{Real}
    vox_vec::Array{Real}
    nb_vox::Vector{Int}
    origin::Vector{Real}
    function Volume()
        new( Array{Real}(1,1,1), Array{Real}(3,3), Vector{Int}(3),Vector{Real}(3) )
    end
    function Volume{T1 <: Real}( nb_vox_iso::T1 )
        if nb_vox_iso > 0
            new( Array{Real}( nb_vox_iso, nb_vox_iso, nb_vox_iso), Array{Real}(3,3), [nb_vox_iso, nb_vox_iso, nb_vox_iso] )
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
    function Volume{T1 <: Real}( nb_vox::Vector{T1} )
        if size(nb_vox)[1] == 3
            new( Array{Real}(nb_vox[1],nb_vox[2],nb_vox[3]),Array{Real}(3,3),nb_vox);
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
    cell_matrix=Array{Real}(3,3);
    for i=1:3
        for j=1:3
            cell_matrix[i,j] = parse(Float64, split( lines[3+i] )[1+j] )*Bohr2Ang
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

    #-----------------------
    # Reads Volumetric Data
    #----------------------------------------------------
    nb_tot=nb_vox[1]*nb_vox[2]*nb_vox[3]
    nb_col=6
    matrix = Array{Real}( nb_vox[1], nb_vox[2], nb_vox[3] )
    offset=(Int)(6+nb_atoms+1)
    x=1; y=1; z=1;
    for i=0:nb_tot/nb_col-1
        for j=1:nb_col
            matrix[x,y,z] = parse(Float64, split( lines[(Int)(offset+i)])[j] )
            z=z+1;
            if z == nb_vox[3]+1
                z=1;
                y=y+1;
            end
            if y == nb_vox[2]+1
                y=1;
                z=1;
                x=x+1;
            end
        end
    end
    #----------------------------------------------------

    #-----------------
    # Updating volume
    #--------------------------------
    volume=Volume()
    volume.nb_vox=nb_vox
    volume.vox_vec=cell_matrix
    volume.matrix=matrix
    volume.origin=center
    #---------------------------------

    #-------------------------------------------
    # Scaling cell vectors by number of voxels
    #--------------------------------------------
    cell_vecs=Cell_matrix(cell_matrix)
    for i=1:3
        cell_vecs.matrix[i,i]=cell_vecs.matrix[i,i]*nb_vox[i]
    end
    #---------------------------------------------

    return atom_list, cell_matrix, volume
end

function getClosest{ T1 <: Real, T2 <: Volume }( position::Vector{T1} , vol::T2 )
    # Works for orthorombic, not yet for non-orthorombic
    #------------------
    # Compute lengths
    #--------------------------------------------
    params=zeros(3,1)
    for i=1:3
        for j=1:3
            params[i] += vol.vox_vec[i,j]
        end
        params[i]=sqrt(params[i])
    end
    # #--------------------------------------------
    #
    # #----------------------------------------------------
    fracs=[0.,0.,0.]
    index=[0,0,0]
    for i=1:3
        print("check: ",(position[i] + vol.origin[i])/params[i], "\n" )
        fracs[i]=( position[i] + vol.origin[i] )/params[i]
        if fracs[i] - trunc(fracs[i]) > 0.5
            index[i]=Int(trunc(fracs[i])+1)
        else
            index[i]=Int(trunc(fracs[i]))
        end
    end
    # #----------------------------------------------------
    #
    # # Returns the index
    return index
end

function paramVoxVectors{ T1 <: Volume }( volume::T1 )
    params=Vector{Real}(3)
    for i=1:3
        for j=1:3
            params[i]=volume.vox_vec[i,j]^2
        end
        params[i] = sqrt(params[i])
    end
    return params
end

# Trace the volume between two points.
function traceLine{ T1 <: Real, T2 <: Real, T3 <: Volume, T4 <: Int }(     position1::Vector{T1}, position2::Vector{T2}, volume::T3, nb_points::T4 )
    if size(position1)[1] != 3 || size(position2)[1] != 3
        return false
    end

    dpos=position2-position1

    # Output Tables
    distances=Vector{Real}(1)

    # Moving along the lines
    curseur=position1
    for i=1:nb_points
        indexs=getClosest(curseur,volume) - position1
        params=paramVoxVectors(volume)
        dist=0
        for j=1:3
            dist += (indexs[i]*params[j])^2
        end
        dist=sqrt(dist)
        value=volume.matrix[indexs[1],indexs[2],indexs[3]]

        curseur += dpos/points
    end

    # Sort by distance from point1

    return dist, values
end

end
