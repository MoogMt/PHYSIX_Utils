function computeCost( distance_matrix::Array{T1,2}, positions::Array{T2}, cost_coeff::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    n_structures=size(distance_matrix)[1]
    cost=0.
    n_dim=size(positions)[2]
    for i=1:n_structures-1
        for j=i+1:n_structures
            dist = 0
            for k=1:n_dim
                dist += ( positions[i,k] - positions[j,k] )*( positions[i,k] - positions[j,k] )
            end
            dist=sqrt(dist)
            cost += 0.5*cost_coeff*( dist - distance_matrix[i,j] )*( dist - distance_matrix[i,j] )
        end
    end
    return cost
end

# Folder where to find file
folder="/home/moogmt/CO2_Classic/"

# Convergence parameters
n_iterations=100000
k=0.02
s=0.005
kT=0.0000002

# Reading file
data=open(string(folder,"FRAME_TO_FRAME.MATRIX"))
number_structure, dmax = split(readline(data))
number_structure=parse(Int64,number_structure)
dmax=parse(Float64,dmax)
distance_matrix=zeros(number_structure,number_structure)
for i=1:number_structure
    print( "Reading FRAME_TO_FRAME: ",i/number_structure*100,"%\n" )
    line=split(readline(data))
    for j=i+1:number_structure
        distance_matrix[i,j]=parse(Float64,line[j])
        distance_matrix[j,i]=distance_matrix[i,j]
    end
end
close(data)

# Start with random positions
x=rand(number_structure)
y=rand(number_structure)

# Compute Initial cost
cost=computeCost(distance_matrix,[x y],k)

for iteration=1:n_iterations
    # Moving randomly a single element
    global element= Int( trunc( rand()*number_structure+1 ) )
    global xx=copy(x)
    global yy=copy(y)
    xx[ element ] = x[element] + ( rand() - 0.5 )*s
    yy[ element ] = y[element] + ( rand() - 0.5 )*s
    # Compute cost of move
    global cost_move=computeCost(distance_matrix,[xx yy],k)
    # Evalutating change of cost
    global r=rand()
    global m=0
    global de=( cost_move - cost )/kT
    if de > 0
        global m=exp(-de)
        if r < m
            global x=copy(xx)
            global y=copy(yy)
            global cost = cost_move
        end
    else
        global x=copy(xx)
        global y=copy(yy)
        global cost = cost_move
    end
end

map_stuff=open(string(folder,"map.dat"),"w")
for i=1:number_structure
    write(map_stuff,string(x[i]," ",y[i],"\n"))
end
close(map_stuff)

error_map=open(string(folder,"error_map.dat"),"w")
for i=1:number_structure
    for j=i+1:number_structure
        write(error_map, string(distance_matrix[i,j]," ",sqrt( ( x[i]-x[j] )*( x[i]-x[j] ) + ( y[i]-y[j] )*( y[i]-y[j] ) ) , "\n") )
    end
end
close(error_map)
