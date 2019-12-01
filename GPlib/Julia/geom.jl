module geom

# Calculate the distance between two vectors
function distance( vec1::Vector{T1}, vec2::Vector{T2} ) where { T1 <: Real , T2 <: Real }
  if size(vec1)[1] != size(vec2)[1]
    error(string("Size mismatch between the two vectors ",size(vec1)[1]," for first argument and",size(vec1)[2],"for second argument.\n"))
  end
  distance=0
  for i=1:size(vec1)[1]
    distance = distance + (vec1[i]-vec2[i])^2
  end
  return sqrt(distance)
end

function norm( vector::Vector{T1} ) where { T1 <: Real }
    norm = 0;
    for i=1:size(vector)[1]
        norm += vector[i]*vector[i];
    end
    return sqrt(norm)
end

function scalar( vector1::Vector{T1} , vector2::Vector{T2} ) where { T1 <: Real, T2 <: Real }
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

function dist2Line( vector_2_line::Vector{T1}, line_vector::Vector{T2} ) where { T1 <: Real, T2 <: Real}
    return abs(scalar(vector_2_line,line_vector))/norm(line_vector)
end

function pointFromLine( point::Vector{T1}, line_vector::Vector{T2} , point_in_line::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    return dist2line( point-point_in_line, line_vector )
end

function angleAlKash( a::T1 , b::T2, c::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    return acosd((a*a+b*b-c*c)/(2*a*b))
end

end
