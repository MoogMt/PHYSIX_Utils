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
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
    function Cell_param()
        new(0,0,0,0,0,0);
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
#---------------------------------------------------------------------------
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
#---------------------------------------------------------------------------

end
