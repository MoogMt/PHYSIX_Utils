include("atoms.jl")

module cell_mod

# Import all import module
#----------------------------
import atom_mod.AtomList
importall atom_mod
#----------------------------

#-------------
# Structures
#-----------------------------
mutable struct Cell_param
    length::Vector{Real}
    angles::Vector{Real}
    function Cell_param()
        new([0,0,0],[0,0,0]);
    end
    function Cell_param{T1 <: Real, T2<: Real, T3 <: Real}( a::T1, b::T2, c::T3)
        new([a,b,c],[90.,90.,90.])
    end
    function Cell_param{ T1 <: Real, T2<: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }( a::T1, b::T2, c::T3, alpha::T4, beta::T5, gamma::T6)
        new([a,b,c],[alpha,beta,gamma])
    end
end
mutable struct Cell_vec
    v1::Vector{Real}
    v2::Vector{Real}
    v3::Vector{Real}
    function Cell_vec()
        new([],[],[]);
    end
end
mutable struct Cell_matrix
    matrix::Array{Real,2}
    function Cell_matrix()
        new(Array{Real}(3,3));
    end
    function Cell_matrix{ T1 <: Real }( matrix::Array{T1})
        if size(matrix)[1]==3 && size(matrix)[2]==3
            new( matrix )
        end
    end
    function Cell_matrix{ T1 <: Real }( a::T1, b::T1, c::T1 )
        new([[a,0,0],[0,b,0],[0,0,c]])
    end
end
#------------------------------

#------------------------------
# General type and conversions
#--------------------------------------------
Cell=Union{Cell_param, Cell_vec, Cell_matrix}
#---------------------------------------------

# Functions
#---------------------------------------------------------------------------\
function vec2matrix{ T1 <: Cell_vec}( vectors::T1 )
    matrix=Cell_matrix()
    for i=1:3
        matrix.matrix[i,1] = vectors.v1[i]
    end
    for i=1:3
        matrix.matrix[i,2] = vectors.v2[i]
    end
    for i=1:3
        matrix.matrix[i,3] = vectors.v3[i]
    end
    return matrix
end
function wrap{ T1 <: Real}( position::T1, length::T1 )
    sign=-1
    if position < 0
        sign=1
    end
    while position < 0 || position > length
        position = position + sign*length
    end
    return position
end
function wrap{ T1 <: AtomList, T2 <: Cell_matrix }( atoms::T1, cell::T2 )
    # Computes cell parameters
    #--------------------------------------------
    params=[0.,0.,0.]
    for i=1:3
        for j=1:3
            params[i]=params[i]+cell.matrix[i,j]
        end
    end
    #--------------------------------------------

    #---------------
    # Compute atoms
    #---------------------------------
    for i=1:size(atoms.positions)[1]
        for j=1:3
             atoms.positions[i,j] = wrap( atoms.positions[i,j],params[j])
        end
    end
    #----------------------------------

    return atoms
end

# Distance related functions
#-------------------------------------------------------------------------------
function dist1D{ T1 <: Real, T2 <: Real, T3 <: Real }( x1::T1, x2::T2, a::T3 )
    dx=x1-x2
    if dx > a*0.5
        return (dx-a)^2
    end
    if dx < -a*0.5
        return (dx+a)^2
    end
    return dx^2
end
function distance{ T1 <: AtomList, T2 <: Cell_param , T3 <: Int }( atoms::T1, cell::T2, index1::T3, index2::T3 )
    distance=0
    for i=1:3
        distance += dist1D( atoms.positions[index1,i],atoms.positions[index2,i], cell.length[i] )
    end
    return sqrt(distance)
end
function distance{ T1 <: AtomList, T2 <: Cell_param,  T3 <: Int, T4 <: Bool }( atoms::T1, cell::T2, index1::T3, index2::T3, wrap::T4 )
    if (  wrap )
        wrap(atoms,cell)
    end
    return distance(atoms,cell,index1,index2)
end
#---------------------------------------------------------------------------

#----------
# Compress
#---------------------------------------------------------------------------
function compressParams{ T1 <: Cell_param, T2 <: Real }( cell::T1, fracs::Vector{T2})
    for i=1:3
        cell.lengths[i] *= fracs[i]
    end
    return cell
end
function compressAtoms{ T1 <: AtomList, T2 <: Cell_param, T3 <: Real }( atoms::T1 , cell::T2, fracs::Vector{T3} )
    for i=1:size(atoms.names)[1]
        for j=1:3
            atoms.positions[i,j] *= fracs[j]
        end
    end
    return atoms
end
#---------------------------------------------------------------------------

end
