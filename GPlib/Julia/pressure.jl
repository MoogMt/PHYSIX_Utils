module pressure

function readPressure( file_name::AbstractString , diag::Bool )
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    nb_lines=size(lines)[1]
    pressure=zeros(nb_lines/4,1)
    if ! diag
        for i=1:nb_lines/4.
            for j=2:4
                pressure[î]= parse(Float64,split(lines[i+j])[j-1])
            end
            pressure[i] = pressure[i]/3.
        end
    else
        pressure_matrix=zeros(3,3)
        for i=1:nb_lines/4
            for j=2:4
                for k=1:3
                    pressure_matrix[j,k]=
                end
            end
        end
    end

    return pressure
end

eigvals(A)

end
