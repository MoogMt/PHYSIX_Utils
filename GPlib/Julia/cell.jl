module cell_mod

mutable struct Cell_param
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
end
mutable struct Cell_vec
    v1::Vector{Real}
    v2::Vector{Real}
    v3::Vector{Real}

print("cell_mod loaded");
end
