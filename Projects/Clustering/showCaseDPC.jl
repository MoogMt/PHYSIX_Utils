GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/Tests/"

function gauss( amplitude::T1, position::Vector{T2}, width::T3,  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    value=0
    for i=1:size(position)[1]
        value += (x[i]-position[i])*(x[i]-position[i])
    end
    return amplitude*exp( - (value)/(2*(width*width)) )
end

function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    value=0
    for i=1:size(amplitudes)[1]
        value += gauss(amplitudes[i],positions[i,:],widths[i],x)
    end
    return value
end

function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Array{T4,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    nb_points=size(x)[1]
    values=zeros(nb_points)
    for i=1:nb_points
        values[i] = gauss( amplitudes, positions,widths,x[i,:])
    end
    return values
end

n_dim=2
seed=trunc(rand()*100)

centers=rand(n_dim,n_dim)
widths=ones(n_dim).*0.1
amplitudes=ones(n_dim)
amplitudes[1]*=3
amplitudes[2]*=1.5


loop_points=[10,50,100,500,1000,2000,5000]
for nb_points in loop_points

    points=rand(nb_points,n_dim)
    value=gauss(amplitudes,centers,widths,points)

    file_out=open(string(folder_base,"data-",nb_points,"-",n_dim,"-",seed,".dat"),"w")
    for i=1:nb_points
        for j=1:n_dim
            write(file_out,string(points[i,j]," "))
        end
        write(file_out,string(value[i],"\n"))
    end
    close(file_out)

    rho_sorted=copy(value)
    index_sorted=clustering.simpleSequence(nb_points)
    for i=1:nb_points
        for j=1:nb_points
            if rho_sorted[i] > rho_sorted[j]
                swap(rho_sorted,i,j)
                swap(index_sorted,i,j)
            end
        end
    end

    delta_sorted=ones(nb_points)*10
    for i=2:nb_points
        for j=1:i-1
            distance=0
            for k=1:n_dim
                distance += (points[index_sorted[i],k]-points[index_sorted[j],k])*(points[index_sorted[i],k]-points[index_sorted[j],k])
            end
            distance=sqrt(distance)
            if delta_sorted[i] > distance
                delta_sorted[i] = distance
            end
        end
    end

    file_out=open(string(folder_base,"decision-diagram-",nb_points,"-",n_dim,"-",seed,".dat"),"w")
    for i=1:nb_points
        write(file_out,string(rho_sorted[i]," ",delta_sorted[i],"\n"))
    end
    close(file_out)

end
