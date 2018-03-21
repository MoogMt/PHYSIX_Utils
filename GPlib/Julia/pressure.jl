module pressure

function readPressureCPMD( file_name::AbstractString , diag::Bool , stride::Int)
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    nb_lines=size(lines)[1]
    nb_pressure_points=Int(trunc(nb_lines/(4*stride)))
    pressure=Vector{Real}(nb_pressure_points)
    for i=1:nb_pressure_points
        pressure[i]=0
    end
    if ! diag
        for i=1:nb_pressure_points
            for j=1:3
                pressure[i] += parse(Float64,split(lines[1+4*(i-1)*stride+j])[j])
            end
            pressure[i] = pressure[i]/3.
        end
    else
        pressure_matrix=zeros(3,3)
        for i=1:nb_pressure_points
            for j=1:3
                line=split(lines[1+4*(i-1)*stride+j])
                for k=1:3
                    pressure_matrix[j,k]=parse(Float64,line[k])
                end
            end
            eigen_press=eigvals(pressure_matrix)
            for l=1:3
                pressure[i] += eigen_press[l]
            end
            pressure[i] = pressure[i]/3.
        end
    end

    return pressure
end

function readStressCPMD( file_name::AbstractString, stride::Int )
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    # Creation de variables
    #----------------------------------------------------------
    nb_stress_points=Int(trunc(size(lines)[1]/(4*stride)))
    stress=Array{Real}(nb_stress_points,6)
    #----------------------------------------------------------
    for i=1:nb_stres_points
        for j=1:3
            for k=1:3
                stress[i,j] = parse(Float64,split(lines[1+4*(i-1)*stride+j])[k])
            end
        end
    end
    #----------------------------------------------------------

    return stress
end

end
