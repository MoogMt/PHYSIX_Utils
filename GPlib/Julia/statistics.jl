module statistics

function simpleAverage{T1 <: Real}( data::Vector{T1})
    average=0.
    size_data=size(data)[1]
    for i=1:size_data
        average += data[i]
    end
    return average/size_data
end

function simpleMoment{T1 <: Real, T2 <: Int}( data::Vector{T1}, n::T2)
    moment=0
    momentn=0
    size_data=size(data)[1]
    for i=1:size_data
        moment += data[i]
        momentn += data[i]^n
    end
    return momentn/size_data - (moment/size_data)^n
end

end
