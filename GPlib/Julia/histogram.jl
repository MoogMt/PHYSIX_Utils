
function minMax( data::Vector{T1} ) where { T1 <: Real }
    min_ = data[1]
    max_ = data[1]
    for i=1:size(data)
        if min_ > data[i]
            min_=data[i]
        end
        if max_ < data[i]
            max_=data[i]
        end
    end
    return min_,max_
end

function histogram( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    min_,max_ = minMax(data)
    histogram = zeros( Int, size(data)[1] )
    for i=1:size(data)[1]
        histogram[ Int(trunc(data[i]/nb_box))+1 ]  += 1
    end
    return histogram
end

function histogramNormed( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    min_,max_ = minMax(data)
    histogram_normed = zeros( Real, size(data)[1] )
    for i=1:size(data)[1]
        histogram[ Int(trunc(data[i]/nb_box))+1 ]  += 1
    end
    histogram /= sum(histogram)
    return histogram
end
