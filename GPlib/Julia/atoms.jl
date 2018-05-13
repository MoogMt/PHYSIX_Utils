include("utils.jl")

module atom_mod

export Atom, AtomList, AtomMolList, switchAtoms

using utils

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
    function AtomList{T1 <: Int}( nb_atoms::T1)
        new( Vector{AbstractString}(nb_atoms) , zeros(nb_atoms), zeros(nb_atoms,3) )
    end
end

mutable struct AtomMolList
    atom_names::Vector{AbstractString}
    atom_index::Vector{Int}
    mol_names::Vector{AbstractString}
    mol_index::Vector{Int}
    positions::Array{Real}
    function AtomMolList{ T1 <: Int }( nb_atoms::T1 )
        new(Vector{AbstractString}( nb_atoms ),Vector{Int}( nb_atoms ) ,Vector{AbstractString}( nb_atoms ), Vector{Int}(nb_atoms), Array{Real}(nb_atoms,3))
    end
end

function switchAtoms{  T2 <: Int, T3 <: Int }( atoms::AtomMolList , index1::T2, index2::T3 )
    # Storing
    a_index=atoms.atom_index[index1]
    a_name=atoms.atom_names[index1]
    m_index=atoms.mol_index[index1]
    m_name=atoms.mol_names[index1]
    positions=atoms.positions[index1,:]
    # Moving 1
    atoms.atom_index[index1]=atoms.atom_index[index2]
    atoms.atom_names[index1]=atoms.atom_names[index2]
    atoms.mol_index[index1]=atoms.mol_index[index2]
    atoms.mol_names[index1]=atoms.mol_names[index2]
    atoms.positions[index1,:]=atoms.positions[index2,:]
    # Moving 2
    atoms.atom_index[index2]=a_index
    atoms.atom_names[index2]=a_name
    atoms.mol_index[index2]=m_index
    atoms.mol_names[index2]=m_name
    atoms.positions[index2,:]=positions
    return
end

end
