module geom

function norm{ T1 <: Real }( vector::Vector{T1} )
    norm = 0;
    for i=1:size(vector)[1]
        norm += vector[i]*vector[i];
    end
    return sqrt(norm)
end

function scalar{ T1 <: Real, T2 <: Real }( vector1::Vector{T1} , vector2::Vector{T2} )
    if size(vector1)[1] != size(vector2)[1]
        # Should return an exception here...
        return false;
    end
    scalar = 0;
    for i=1:size(vector1)[1]
        scalar += vector1[i]*vector2[i]
    end
    return scalar
end

function dist2Line{ T1 <: Real, T2 <: Real}( vector_2_line::Vector{T1}, line_vector::Vector{T2} )
    return abs(scalar(vector_2_line,line_vector))/norm(line_vector)
end

function pointFromLine{ T1 <: Real, T2 <: Real, T3 <: Real }( point::vector{T1}, line_vector::Vector{T2} , point_in_line::Vector{T3} )
    return dist2line( point-point_in_line, line_vector )
end

end
