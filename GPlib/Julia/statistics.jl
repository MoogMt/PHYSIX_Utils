module statistics

function simpleAverage{T1 <: Real}( data::Vector{T1})
    average=0
    size=size(data)[1]
    for i=1:size
        average += data[i]
    end
    return average /= size
end

function simpleMoment{T1 <: Real, T2 <: Int}( data::Vector{T1} , n::T2)
    moment=0
    size=size(data)[1]
    for i=1:size
        moment += data[i]
        momentn += data[i]^n
    end
    return (moment^n - momentn)/size
end

end
