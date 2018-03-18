module pressure

function readPressure( file_name::AbstractString , diag::Bool , stride::Int)
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    nb_lines=size(lines)[1]
    nb_pressure_points=nb_lines/(4*stride)
    pressure=zeros(nb_pressure_points,1)
    if ! diag
        for i=1:nb_pressure_points
            for j=2:4
                pressure[î]= parse(Float64,split(lines[i*stride+j])[j-1])
            end
            pressure[i] = pressure[i]/3.
        end
    else
        pressure_matrix=zeros(3,3)
        for i=1:nb_pressure_points
            for j=2:4
                for k=1:3
                    pressure_matrix[j,k]=split(lines[i*stride+j])[k]
                end
            end
            eigen_press=eigvals(pressure_matrix)
            for l=1:3
                pressure[i]=eigen_press[l]
            end
            pressure[i] = pressure[i]/3.
        end
    end

    return pressure
end

eigvals(A)

end
