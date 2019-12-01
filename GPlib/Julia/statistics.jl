module statistics

using Statistics
using Bootstrap

function average( data::Vector{T1} ) where { T1 <: Real }
    average=0
    nb_point=size(data)[1]
    for i=1:nb_point
        average += data[i]
    end
    return average/nb_point
end

function variance( data::Vector{T1} ) where { T1 <: Real }
    variance=0
    average=0
    nb_point=size(data)[1]
    for i=1:nb_point
        average += data[i]
        variance += data[i]*data[i]
    end
    average /= nb_point
    return variance/nb_point - average*average
end

function averageVariance( data::Vector{T1} ) where { T1 <: Real }
    variance=0
    average=0
    nb_point=size(data)[1]
    for i=1:nb_point
        average += data[i]
        variance += data[i]*data[i]
    end
    average /= nb_point
    return average, variance/nb_point - average*average
end

function blockAverage( data::Vector{T1}, size_block::T2 ) where { T1 <: Real, T2 <: Int }
    
end

function blockAverage( data::Vector{T1}, min_block_size::T2, max_block_size::T3, stride_block::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int }
    nb_point=size(data)[1]

    if min_block_size <= 1
        return zeros(nb_point),zeros(nb_point),zeros(nb_point)
    end

    blocks_number=Int((max_block_size-min_block_size)/stride_block)+1

    meanBlocks=zeros(blocks_number)
    varBlocks=zeros(blocks_number)
    blocks=zeros(blocks_number)

    block_ctrl=1
    for block_size=min_block_size:stride_block:max_block_size
        n_block=Int(floor(nb_point/block_size))
        observation_block=zeros(n_block)
        for i=1:n_block
            beg=(i-1)*block_size+1
            observation_block[i]=Statistics.mean(data[beg:beg+block_size-1])
        end
        meanBlocks[block_ctrl] = Statistics.mean(observation_block)
        varBlocks[block_ctrl] = Statistics.var(observation_block)/(n_block-1)
        blocks[block_ctrl]=block_size
        block_ctrl += 1
    end
    return blocks, meanBlocks, varBlocks
end

function bootStrap( data::Vector{T1}, n_boot::T2) where { T1 <: Real, T2 <: Int }
    bs1 = bootstrap(Statistics.mean, data, Bootstrap.BalancedSampling(n_boot))
    return bs1.t0[1], stderror(bs1)[1]
end

function minMax( data::Vector{T1} ) where { T1 <: Real }
    return minimum(data), maximum(data)
end

function minMax( data::Array{T1,2} ) where { T1 <: Real }
    nb_=size(data)[2]
    mins_=zeros(nb_)
    maxs_=zeros(nb_)
    for dim=1:nb_
        mins_[dim] = minimum( data[:,dim] )
        maxs_[dim] = maximum( data[:,dim] )
    end
    return mins_,maxs_
end

function histogram( data::Vector{T1}, nb_box::T2, min_::T3, max_::T4 ) where { T1 <: Real , T2 <: Int, T3 <: Real, T4 <: Real }
    histogram = zeros( nb_box )
    delta_box=(max_-min_)/nb_box
    for i=1:size(data)[1]
        histogram[ Int(trunc((data[i]-min_)/delta_box))+1 ]  += 1
    end
    return histogram
end

function histogram( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    min_,max_ = minMax(data)
    return histogram( data, nb_box, min_, max_ )
end

function histogramNormed( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    histogram = histogram( data, nb_box )
    histogram /= sum(histogram)
    return histogram
end

function histogramNormed( data::Vector{T1}, nb_box::T2, min_::T3, max_::T4 ) where { T1 <: Real , T2 <: Int, T3 <: Real, T4 <: Real }
    histogram = histogram( data, nb_box, min_, max_ )
    histogram /= sum(histogram)
    return histogram
end

function writeHistogram( file_out::T1, histogram::Vector{T2}, nb_box::T3, min_::T4, max::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Real }
    delta_box=(max_-min_)/nb_box
    file_o = open(file_out, "w")
    for box=1:nb_box
        Base.write(file_o,string(box*delta_box+min_," ",histogram[box],"\n"))
    end
    close(file_o)
    return
end

function histogram2D( data::Array{T1,2}, nb_box::Vector{T2}, mins_::Vector{T3}, maxs::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    delta_box = zeros(Int,2)
    for dim=1:2
        delta_box[dim] = Int( trunc( (maxs_[dim]-mins_[dim])/nb_box[dim] )  )
    end
    histogram = zeros( nb_box[1], nb_box[2] )
    for point=1:nb_point
        histogram[ Int(trunc( (data[point,1]-mins_[1])/delta_box[1] ))+1 , Int(trunc( (data[point,2]-mins_[2])/delta_box[2] ))+1 ]  += 1
    end
    return histogram
end

function histogram2D( data::Array{T1,2}, nb_box::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    mins_,maxs_ = minMax(data)
    return histogram2D( data, nb_box, mins_, maxs_)
end

function histogram2DNormed( data::Array{T1,2}, nb_box::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    histogram = histogram2D( data, nb_box )
    histogram /= sum(histogram)
    return histogram
end

function histogram2DNormed( data::Array{T1,2}, nb_box::Vector{T2}, mins_::Vector{T3}, maxs_::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    histogram = histogram2D( data, nb_box, mins_, maxs_)
    histogram /= sum(histogram)
    return histogram
end

function writeHistogram2D( file_out::T1, histogram::Vector{T2}, nb_box::T3, mins_::T4, maxs_::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Real }
    delta_box=zeros(2)
    for i=1:2
        delta_box[i] = (maxs_[i]-mins_[i])/nb_box[i]
    end
    file_o = open(file_out, "w")
    for box1=1:nb_box[1]
        for box2=1:nb_box[2]
            Base.write(file_o,string(box1*delta_box[1]+mins_[1]," "))
            Base.write(file_o,string(box2*delta_box[2]+mins_[2]," "))
            Base.write(file_o,string(histogram[box1,box2],"\n"))
        end
    end
    close(file_o)
    return
end

end
