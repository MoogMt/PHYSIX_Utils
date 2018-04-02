module atom_mod

mutable struct Atom
    name::AbstractString
    index::Int
    position::Vector{Real}
    function Atom()
        new("",-1,[0,0,0]);
    end
end

mutable struct AtomList
    names::Vector{AbstractString}
    index::Vector{Int}
    positions::Array{Real}
    function AtomList()
        new(Vector{Real}(),Vector{Real}(),Array{Real}(0,1))
    end
    function AtomList{T1 <: Int}( nb_atoms::T1)
        new( Vector{AbstractString}(nb_atoms), Vector{Int}(nb_atoms), Array{Real}(nb_atoms,3) )
    end
end

mutable struct Mass
    mass::Real
end

mutable struct MassList
    mass::Vector{Real}
end

mutable struct Charge
    charge::Real
end

mutable struct ChargeList
    charge::Vector{Real}
end

mutable struct Pp
    pp_name::AbstractString
    pp_path::AbstractString
end

mutable struct PPList
    pp_name::Vector{AbstractString}
    pp_path::Vector{AbstractString}
end

mutable struct Molecule
    name::AbstractString
    index::Int
    atoms::AtomList
end

mutable struct MoleculeList
    name::Vector{AbstractString}
    index::Vector{Int}
    atoms::Vector{AtomList}
end

function distance{ T1 <: Real, T2 <: Real }( atom1::Vector{Real}, atom2::Vector{Real} )
    dist=0
    for i=1:3
        dist+=(atom1[i]-atom2[i])^2
    end
    return sqrt(dist)
end

function distance{ T1 <: Atom, T2 <: Atom }( atom1::T1, atom2::T2 )
    dist=0
    for i=1:3
        dist += (atom1.position[i]-atom2.position[i])^2
    end
    return sqrt(dist)
end

function distance{ T1 <: AtomList, T2 <: Int, T3 <: Int}( atoms::T1, index1::T2, index2::T3 )
    dist=0
    for i=1:3
        dist += ( atoms.positions[index1,i] - atoms.positions[index2,i])^2
    end
    return sqrt(dist)
end

end
