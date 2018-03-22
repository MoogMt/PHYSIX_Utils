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
    function Histogram{ T1 <: Real, T2 <: Int }( data::Vector{T1}, nb_box::T2 )
        # Min and Max
        min=min(data)
        max=max(data)
        # Size of boxes
        delta=(max-min)/delta
        # Construction of the boxes
        begins=Vector{Real}(nb_box)
        ends=Vector{Real}(nb_box)
        for i=1:nb_box
            begins[i]=min
            ends[i]=
            min +=delta
        end
        new(begins,ends,value)
    end
end
#-----------------------------------------------------------------------

end
