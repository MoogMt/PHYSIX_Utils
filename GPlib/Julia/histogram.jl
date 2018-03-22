module histogram

#-----------------------------------------------------------------------
mutable struct Bin
    begin_box::Real
    end_box::Real
    value::Real
end

mutable struct Histogram
    begins_box::Vector{Real}
    ends_box::Vector{Real}
    values::Vector{Real}
    function Histogram()
        new( [],[],[] )
    end
    function Histogram{ T1 <: Real, T2 <: Int }( data::Array{T1}, nb_box::T2 )
        if size(data)[2] == 1
            
        else
            print("Error on size of data array. \n")
            return
        end
    end
end
#-----------------------------------------------------------------------

end
