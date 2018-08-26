# Loading file
include("contactmatrix.jl")

# Folder where to find file
folder="/home/moogmt/CO2/CO2_AIMD/clusters/"

# Convergence parameters
n_iterations=5000000
k=0.02
s=0.005
kT=0.0000002

# Others
best_max_distance = 0.

# Reading frame to frame matrix
if isfile( string(folder,"FRAME_TO_FRAME.MATRIX") )

    # Reading file
    data=open(string(folder,"FRAME_TO_FRAME.MATRIX"))
    number_structure, dmax = split(readline(data))
    number_structure=parse(Int64,number_structure)
    dmax=parse(Float64,dmax)
    d_real=zeros(number_structure,number_structure)
    for i=1:number_structure
        print( "Reading FRAME_TO_FRAME: ",i/number_structure*100,"%\n" )
        line=split(readline(data))
        for j=i+1:number_structure
            d_real[i,j]=parse(Float64,line[j])
            d_real[j,i]=d_real[i,j]
        end
    end
    close(data)

    # Start with random positions
    x=rand(number_structure)
    y=rand(number_structure)

    # Compute and print Initial energy
    energy=0
    for i=1:number_structure-1
        for j=i+1:number_structure
            dij=sqrt( (x[i]-x[j] )^2 + ( y[i]-y[j])^2 )
            energy += 0.5*k*( dij - d_real[i,j] )^2
        end
    end
    print("Initial energy: ",energy,"\n")

    # Copy of position
    xx=zeros(number_structure)
    yy=zeros(number_structure)
    for i=1:number_structure
        xx[i]=x[i]
        yy[i]=y[i]
    end

    for iteration=1:n_iterations
        print("Progress: ",iteration/n_iterations*100,"%\n")
        # Moving randomly a single element
        element= Int( trunc( rand()*number_structure+1 ) )
        xx[ element ] = x[element] + ( rand() - 0.5 )*s
        yy[ element ] = y[element] + ( rand() - 0.5 )*s
        energy_bis=0
        max_distance=0
        for i=1:number_structure-1
            for j=i+1:number_structure
                d_map=sqrt( ( xx[i]-xx[j] )^2 + ( yy[i]-yy[j] )^2 )
                dd = ( d_map - d_real[i,j] )^2
                if dd > max_distance
                    max_distance = dd
                end
                energy_bis += 0.5*k*dd
                if iteration == n_iterations
                    print("error between ",i," and ",j," : ",(d_map-d_real[i,j]), ", relative: ",(d_map-d_real[i,j])/d_real[i,j],", orig_d mapped_d: ",d_map," ",d_real[i,j],"\n")
                end
            end
        end
        # Evalutating change of energy
        r=rand()
        m=0
        de=( energy_bis - energy )/kT
        if de > 0
            m=exp(-de)
        else
            m = r+1
        end
        if r < m
            # Metropolis choice
            for i=1:number_structure
                x[i] = xx[i]
                y[i] = yy[i]
            end
            energy = energy_bis
        end
    end

    map=open(string(folder,"map.dat"),"w")
    for i=1:number_structure
        write(map,string(x[i]," ",y[i],"\n"))
    end
    close(map)

    error_map=open(string(folder,"error_map.dat"),"w")
    for i=1:number_structure
        for j=i+1:number_structure
            write(error_map, string(d_real[i,j]," ",sqrt( ( x[i]-x[j] )^2 + ( y[i]-y[j] )^2 ) , "\n") )
        end
    end
    close(error_map)
else
    print("OUPS")
end
