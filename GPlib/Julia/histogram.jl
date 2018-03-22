module histogram

#-----------------------------------------------------------------------
mutable struct bin
    begin_box::Real
    end_box::Real
    value::Real
end

mutable struct histogram
    begins_box::Vector{Real}
    ends_box::Vector{Real}
    values::Vector{Real}
    function Volume()
        new( Array{Real}(0,4) )
    end
end
#-----------------------------------------------------------------------

end
