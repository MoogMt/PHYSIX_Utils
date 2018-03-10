module cell_mod

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
    matrix::Array{Real}
    function Cell_matrix()
        new([]);
    end
end

print("cell_mod loaded");

end
