module cell_mod

include("atoms.jl")

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
mutable struct Cell_vec
    v1::Vector{Real}
    v2::Vector{Real}
    v3::Vector{Real}
    function Cell_vec()
        new([],[],[]);
    end
end
end
mutable struct Cell_matrix
    matrix::Array{Real,2}
    function Cell_matrix()
        new(Array{Real}(3,3));
    end
    function Cell_matrix{ T1 <: Real }( a::T1, b::T1, c::T1 )
        new([[a,0,0],[0,b,0],[0,0,c]])
    end
end

Cell=Union{Cell_param, Cell_vec, Cell_matrix}

function wrap{ T1 <: Real}( position::T1, length::T1 )
    sign=1
    if position < 0
        sign=-1
    end
    while position < 0 || position > length
        position = position + sign*length
    end
    return position
end

function wrap{ T1 <: atom_mod.AtomList, T2 <: Cell_matrix }( atoms::T1, cell::T2 )
    # Computes cell parameters
    #--------------------------------------------
    param=[0.,0.,0.]
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
            atoms.position[i,j] = wrap( atoms.positions[i,j],param[i])
        end
    end
    #----------------------------------

    return atoms
end

print("cell_mod loaded");

end
