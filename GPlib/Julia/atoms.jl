include("utils.jl")

print("Loading Atoms")

module atom_mod

mutable struct Atom
    name::AbstractString
    index::Int
    position::Vector{Real}
    function Atom()
        new("",-1,[0,0,0]);
    end
end
export Atom

mutable struct AtomList
    names::Vector{AbstractString}
    index::Vector{Int}
    positions::Array{Real}
    function AtomList{T1 <: Int}( nb_atoms::T1)
        new( Vector{AbstractString}(nb_atoms) , zeros(nb_atoms), zeros(nb_atoms,3) )
    end
end
export AtomList

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
export AtomMolList

function distance{ T1 <: Real, T2 <: Real }( atom1::Vector{T1}, atom2::Vector{T2} )
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
export distance

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
export switchAtoms

end
