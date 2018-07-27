# Loading file
include("contactmatrix.jl")

# Folder where to find file
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"

# number of centers
number_structure = 20

# Convergence parameters
n_iterations=100000000
k=0.05
s=0.002
kT=0.00000001

# Others
best_max_distance = 0.

# Reading frame to frame matrix
if isfile( string(folder,"FRAME_TO_FRAME.MATRIX") )

    # Reading file
    data=open(string(folder,"FRAME_TO_FRAME.MATRIX"))
    map_check=open(string(folder,"map_check"),"w")
    n, dmax = split(readline(data))
    n=parse(Int64,n)
    dmax=parse(Float64,dmax)
    d_real=zeros(n,n)
    for i=1:n
        print( "Reading FRAME_TO_FRAME: ",i/n*100,"%\n" )
        line=split(readline(data))
        for j=1:n
            d_real[i,j]=parse(Float64,line[j])*dmax
            write(map_check,string(i," ",j," ",d_real[i,j],"\n"))
        end
    end
    close(data)
    close(map_check)

    # Start with random positions
    pos=rand(number_structure,2)

    # Compute and print Initial energy
    energy=0
    for i=1:number_structure-1
        for j=i+1:number_structure
            dij=sqrt( (pos[i,1]-pos[j,1])^2 + (pos[i,2]-pos[j,2])^2 )
            energy += 0.5*k*( dij - d_real[i,j] )^2
        end
    end

    # Copy of position
    pos_bis=pos

    error_map=open(string(folder,"error_map.dat"),"w")
    conv=open(string(folder,"conv.dat"),"w")
    for iteration=1:n_iterations
        print("Progress: ",iteration/n_iterations*100,"%\n")
        # Moving randomly a single element
        element= Int( trunc( rand()*number_structure+1 ) )
        pos_bis[ element , 1 ] = pos_bis[ element , 1 ] + ( rand() - 0.5 )*s
        pos_bis[ element , 2 ] = pos_bis[ element , 2 ] + ( rand() - 0.5 )*s
        energy_bis=0
        max_distance=0
        for i=1:number_structure-1
            for j=i+1:number_structure
                d_map=sqrt( ( pos_bis[i,1]-pos_bis[j,1] )^2 + ( pos_bis[i,2]-pos_bis[j,2] )^2 )
                dd = sqrt( ( d_map - d_real[i,j] )^2 )
                if dd > max_distance
                    max_distance = d_map
                end
                energy_bis += 0.5*k*dd
                if iteration == n_iterations
                    print("error between ",i," and ",j," : ",(d_map-d_real[i,j]), ", relative: ",(d_map-d_real[i,j])/d_real[i,j],", orig_d mapped_d: ",d_map," ",d_real[i,j],"\n")
                    write(error_map,string(d_real[i,j]," ",d_map,"\n"))
                end
            end
        end
        # Evalutating change of energy
        r=rand()
        de=( energy_bis - energy )/kT
        if de > 0
            if de < 50
                m=exp(-de)
            else
                m=0
            end
        else
            m=2
        end
        if r < m
            # Metropolis choice
            for i=1:number_structure
                pos[i,:] = pos_bis[i,:]
            end
            energy = energy_bis
        end
        best_max_distance=max_distance
        write(conv,string(iteration," ",max_distance," ",energy,"\n"))
    end
    close(error_map)
    close(conv)

    print("Final energy: ",energy," average and max deviation:",sqrt(2.*e/(k*n*(n-1)/2.))," ",sqrt(best_max_distance),"\n" )

    map=open(string(folder,"map.dat"),"w")
    for i=1:number_structure
        for j=1:2
            write(map,string(pos[i,j]," "))
        end
        write(map,"\n")
    end
    close(map)
else
    print("OUPS")
end
