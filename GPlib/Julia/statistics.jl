module statistics

function simpleAverage( data::Vector{T1}) where {T1 <: Real}
    average=0.
    size_data=size(data)[1]
    for i=1:size_data
        average += data[i]
    end
    return average/size_data
end

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

function blockAverage( data::Vector{T1}, block_size::T2 ) where {T1 <: Real,  T2 <: Int}
    average=0
    data_size=size(data)[1]
    nb_block=Int(trunc(data_size/block_size))
    # Compute the average of averages
    count=0
    for i=1:nb_block
        if data_size > block_size+i*block_size
            # Compute the average inside the group
            local_average=0
            for j=1:block_size
                local_average += data[i*block_size+j]
            end
            average += local_average/block_size
            count = count + 1
        end
    end
    return average/count
end

function blockMoment( data::Vector{T1}, block_size::T2 , n::T3) where { T1 <: Real, T2 <: Int, T3 <: Int }
    moment=0; momentn=0;
    data_size = size(data)[1]
    nb_block = Int(trunc(size_data/block_size))
    # Compute the variance of variance
    for i=1:nb_block
        local_moment=0; local_moment=0;
        for j=1:block_size
            local_moment += local_moment
            local_momentn += local_momentn^n
        end
        variance=local_momentn/block_size - (local_moment/block_size)^n
        moment += variance
        momentn += variance^n
    end
    return momentn/nb_block - (moment/nb_block)^n
end

function bootStrap( data::Vector{T1}, sample_size::T2, n::T2) where { T1 <: Real, T2 <: Int}
    average=0; moment=0;
    return average, moment;
end

end
