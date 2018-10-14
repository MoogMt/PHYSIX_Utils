include("contactmatrix.jl")

function gaussian{ T1 <: Real, T2 <: Real, T3 <: Real}( position::T1, center::T2, width::T3)
    return exp( -(center-position)*(center-position)/(2*width*width))
end

function gaussian{ T1 <: Real, T2 <: Real, T3 <: Real }( position::Vector{T1}, center::Vector{T2}, width::Vector{T3} )
    value=0
    for i=1:size(position)[1]
        value += gaussian(position[i],center[i],width[i])
    end
    return value
end
