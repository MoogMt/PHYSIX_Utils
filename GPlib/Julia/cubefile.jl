include("cell.jl")
include("geom.jl")

module cube_mod

include("conversion.jl");

import atom_mod.AtomList
import cell_mod.Cell_matrix
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
    function Volume()
        new( Array{Real}(0,4) )
    end
    function Volume{T1 <: Real}( nb_vox_iso::T1 )
        if nb_vox_iso > 0
            new( Array{Real}(nb_vox_iso*nb_vox_iso*nb_vox_iso,4),Array{Real}(3,3),[nb_vox_iso,nb_vox_iso,nb_vox_iso])
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
    function Volume{T1 <: Real}( nb_vox::Vector{T1} )
        if size(nb_vox)[1] == 3
            nb_tot = nb_vox[1]*nb_vox[2]*nb_vox[3];
            new( Array{Real}(nb_tot,4),nb_vox);
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
            cell_matrix.matrix[i,j] = parse(Float64, split( lines[3+i] )[1+j] )*Bohr2Ang
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
    volume = Volume(nb_vox)
    volume.nb_vox=nb_vox
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

function getClosest{ T1 <: Real}( position::Vector{T1} , volume::Volume )
    # Works for orthorombic
    params=[0,0,0]
    # Compute lengths
    for i=1:3
        for j=1:3
            params[i]=volume.vox_vec[i,j]**2
        end
        params[i]=sqrt(params[i])
    end
    indexs=[0,0,0]
    for i=1:3
        indexs[i]=position[i]/params[i]
        if indexs[i] - trunc(indexs[i]) > 0.5
            indexs[i]=
        end
    end
    # Returns the index
    return [position[1]/params[1],position[2]/params[2],position[3]/params[3]]
end

# Trace the volume between two points.
function traceVolume{ T1 <: Real, T2 <: Real, T3 <: Volume }( position1::Vector{T1}, position2::Vector{T2}, volume::T3 )
    if size(position1)[1] != 3 || size(position2)[1] != 3
        return false
    end
    indexs1= getClosest( position1, volume);
    indexs2= getClosest( position2, volume);
    # dindex=index1-index2
    # Getting the possible direction
    # direction=Vector{Int}(3)
    # for i=1:3
    #     if dindex[i] > 0
    #         direction[i] = 1
    #     elseif dindex[i] < 0
    #         direction[i] = -1
    #     else
    #         direction = 0
    #     end
    # end
    # # Making the move matrix
    # moveMatrix=Array{Int}(7,3)
    # for i=0:1
    #     for j=0:1
    #         for k=0:1
    #             if i+j+k != 0
    #                 moveMatrix[i,:]=[i*direction[1],j*direction[2],k*direction[3]])
    #             end
    #         end
    #     end
    # end
    # distances=Array{Real}(7,1)
    # Clearing direction
    # curseur=indexs1
    list=Array{Real}(0,3)
    # while norm(curseur-indexs2) > 0
    #     for i=1:7
    #         test=curseur+moveMatrix[i,:]
    #         distance2line(test,)
    #     end
    #     vcat(index2,curseur)
    # end
end

print("Cube Module Loaded!\n")

end
