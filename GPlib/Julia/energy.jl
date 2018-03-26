module energy

function computeH{ T1 <: Real, T2 <: Real, T3 <: Real }( energy::Vector{T1},  pressure::Vector{T2} ,volume::Vector{T3})
    sU=size(energy)[1]
    sP=size(energy)[1]
    sV=size(energy)[1]
    if sU == sP && sU == sV
        
    else
        print("Error: Vectors do not have the same sizes: E->",sU," P->",sP," V->",sU)
        return false
    end
end

end
