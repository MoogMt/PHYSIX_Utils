module contact_matrix

mutable struct ContactMatrix
    ContactMatrix::Array{Real,2}
end

mutable struct ContactVector
    ContactVector::Vector{Real}
end

end
