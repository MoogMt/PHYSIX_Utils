module statistics

using Statistics
using Bootstrap

function simpleMoment( data::Vector{T1}, n::T2) where {T1 <: Real, T2 <: Int}
    moment=0
    momentn=0
    size_data=size(data)[1]
    for i=1:size_data
        moment += data[i]
        momentn += data[i]^n
    end
    return momentn/size_data - (moment/size_data)^n
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

end
