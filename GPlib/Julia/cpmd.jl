module cpmd

using conversion

export readInputTimestep
export readEnergy, readPressure, readStress

function readInputTimestep( file_input_path::T1 ) where { T1 <: AbstractString }
    # Read input
    file_in=open(file_input_path)
    lines=readlines(file_in)
    close(file_in)

    # Extract timestep
    timestep=0
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "TIMESTEP"
                timestep=parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end

    # Conversion to fs
    return conversion.Hatime2fs*timestep
end

function strideStress( file_input_path::T1 ) where { T1 <: AbstractString }
    # Readinput
    file_in=open(file_input_path)
    lines=readlines(file_in)
    close(file_in)

    # Extract stride of STRESS tensor
    stride_stress = 0
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "STRESS TENSOR"
                stride_stress = parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end

    return stride_stress
end

function readEnergy( file_name::T1 ) where { T1 <: AbstractString }
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    # Array Init
    #----------------------------------------
    nb_steps=size(lines)[1]
    temperature=Vector{Real}(undef,nb_steps)
    e_class=Vector{Real}(undef,nb_steps)
    e_ks=Vector{Real}(undef,nb_steps)
    msd=Vector{Real}(undef,nb_steps)
    time=Vector{Real}(undef,nb_steps)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    for i=1:nb_steps
        line=split(lines[i])
        temperature[i]=parse(Float64,line[3])
        e_ks[i]=parse(Float64,line[4])
        e_class[i]=parse(Float64,line[5])
        msd[i]=parse(Float64,line[7])
        time[i]=parse(Float64,line[8])
    end
    #----------------------------------------------

    return  temperature, e_ks, e_class, msd, time
end

function readPressure( file_name::T1 , diag::T2 , stride::T3 ) where { T1 <: AbstractString, T2 <: Bool, T3 <: Int }
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    nb_lines=size(lines)[1]
    nb_pressure_points=Int(trunc(nb_lines/(4*stride)))
    pressure=Vector{Real}(undef,nb_pressure_points)
    for i=1:nb_pressure_points
        pressure[i]=0
    end
    if ! diag
        for i=1:nb_pressure_points
            for j=1:3
                if split(lines[1+4*(i-1)*stride+j])[j] == "TOTAL"
                    print("Target is line ", 1+4*(i-1)*stride+j)
                    print(" CHECK TOTAL\n")
                    print("element=",j,"\n")
                    print(lines[1+4*(i-1)*stride+j],"\n")
                    quit()
                end
                if split(lines[1+4*(i-1)*stride+j])[j] == "STRESS"
                    print("Target is line ", 1+4*(i-1)*stride+j)
                    print(" CHECK STRESS\n")
                    print("element=",j,"\n")
                    print(lines[1+4*(i-1)*stride+j],"\n")
                    quit()
                end
                if split(lines[1+4*(i-1)*stride+j])[j] == "STEP:"
                    print("Target is line ", 1+4*(i-1)*stride+j)
                    print(" CHECK STEP\n")
                    print("element=",j,"\n")
                    print(lines[1+4*(i-1)*stride+j],"\n")
                    quit()
                end
                if split(lines[1+4*(i-1)*stride+j])[j] == "TENSOR"
                    print("Target is line ", 1+4*(i-1)*stride+j)
                    print(" CHECK TENSOR\n")
                    print("element=",j,"\n")
                    print(lines[1+4*(i-1)*stride+j],"\n")
                    quit()
                end
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

function readStress( file_name::T1, stride::T2 ) where { T1 <: AbstractString, T2 <: Int }
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
