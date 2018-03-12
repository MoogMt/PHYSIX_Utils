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
        new([],[],[])
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

end
