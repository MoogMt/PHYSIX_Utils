#==========#
# Laio Algo
#==============================================================================#
include("clustering.jl")
include("contactmatrix.jl")

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
stride = 1
unit=ps2fs*timestep*stride
start_time=5
start_step=Int(start_time/unit)
nbC=32
nbO=2*nbC

#folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
folder=string("/home/moogmt/CO2/CO2_AIMD/8.82/3000K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.read( string(folder,file), stride, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
maxN=4
n_dim=maxN+1
data_set=zeros(nb_steps*nbC,n_dim)
fileC=open(string(folder,"distancesNN.dat"),"w")
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
		for i=1:maxN
			data_set[ carbon+nbC*(step-1), i ] = distances[ i ]
            write(fileC, string(distances[i]," ") )
		end
        a=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[1]))
        b=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[2]))
        c=cell_mod.distance(traj[step],cell,nbC+Int(index[1]),nbC+Int(index[2]))
        angle=acosd((a*a+b*b-c*c)/(2*a*b))
        data_set[ carbon+nbC*(step-1), maxN+1 ] = angle
        write(fileC, string(angle," ") )
        write(fileC,string("\n"))
	end
end
close(fileC)

size_data=8000
data_train=data_set[1:size_data,:]

max=data_train[1,:]
min=data_train[1,:]

for i=1:size_data
    for j=1:n_dim
        if max[j] < data_train[i,j]
            max[j] = data_train[i,j]
        end
        if min[j] > data_train[i,j]
            min[j] = data_train[i,j]
        end
    end
end


print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim, max, min )

using PyPlot

dc=0.1

function gaussianKernel{ T1 <: Real, T2 <: Real }( matrix_distance::Array{T1,2} , cut_off_distance::T2)
    nb_element=size(matrix_distance)[1]
    rho=zeros(nb_element)
    for i=1:nb_element
        for j=i+1:nb_element
            gauss = exp(-(matrix_distance[i,j]/cut_off_distance)*(matrix_distance[i,j]/cut_off_distance))
            rho[i] += gauss
            rho[j] += gauss
        end
    end
    return rho
end

rho=gaussianKernel(distance_matrix,dc)

function max_array{ T1 <: Real }( distance_matrix::Array{T1,2} )
    nb_element=size(distance_matrix)[1]
    max=0
    for i=1:nb_element
        for j=i+1:nb_element
            if distance_matrix[i,j] > max
                max=distance_matrix[i,j]
            end
        end
    end
    return max
end

max_distance=max_array(distance_matrix)

function simpleSequence{ T1 <: Int }( size_::T1 )
    vector=zeros(Int, size_)
    for i=1:size_
        vector[i] = i
    end
    return vector
end

function sort_descend{ T1 <: Real }( vector::Vector{T1} )
    size_vector=size(vector)[1]
    vector_sorted=copy(vector)
    index_vector=simpleSequence(size_vector)
    for i=1:size_vector
        for j=i+1:size_vector
            if vector_sorted[i] < vector_sorted[j]
                stock=vector_sorted[i]
                stock2=index_vector[i]
                vector[i]=vector_sorted[j]
                index_vector[i]=index_vector[j]
                vector_sorted[j]=stock
                index_vector[j]=stock2
            end
        end
    end
    return vector_sorted, index_vector
end

rho_sorted, index_rho = sort_descend(rho)
# rho=zeros(size_data)
# for i=1:size_data
#     for j=1:size_data
#         if distance_matrix[i,j] < dc && i != j
#             rho[i] += 1
#         end
#     end
# end

delta=zeros(size_data)
delta[ index_rho[1] ] = -1
nneigh=zeros(Int,size_data)
for i=2:size_data
    delta[ index_rho[i] ] = max_distance
    for j=1:i-1
        if distance_matrix[ index_rho[i] , index_rho[j] ] < delta[ index_rho[i] ]
            delta[ index_rho[i] ] = distance_matrix[ index_rho[i], index_rho[j] ]
            nneigh[ index_rho[i] ] = index_rho[j]
        end
    end
end

function max_vector{ T1 <: Real }( vector::Vector{T1} )
    max=0
    for i=1:size(vector)[1]
        if max < vector[i]
            max=vector[i]
        end
    end
    return max
end

delta[ index_rho[1] ] = max_vector( delta )

using PyPlot
PyPlot.figure(1)
PyPlot.plot(rho,delta,"r.")
PyPlot.xlabel("rho")
PyPlot.ylabel("delta")
PyPlot.show()

min_delta=0.119
min_rho=20

n_cluster=0
cl=ones(Int,size_data)
cl*=-1
icl=[]
for i=1:size_data
    if rho[i] > min_rho && delta[i] > min_delta
        n_cluster+=1
        cl[i] = n_cluster
        icl=push!(icl,i)
    end
end

for i=1:size_data
    if cl[ index_rho[i] ] == -1
        cl[ index_rho[i] ] = cl[ nneigh[ index_rho[i] ] ]
    end
end

for i=1:n_cluster
    file=open(string(folder,"cluster-",i,"-index.dat"),"w")
    for j=1:size_data
        if i == cl[j]
            for k=1:n_dim
                write(file,string(data_train[ j, k ]," ",))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end


# Do not know exacly what this is...
halo=copy(cl)

bord_rho=zeros(size_data)
for i=1:size_data-1
    for j=i+1:size_data
        if ( cl[i] == cl[j] ) && ( distance_matrix[i,j] <= dc )
            rho_average = (rho[i]+rho[j])/2
            if rho_average > bord_rho[ cl[i] ]
                bord_rho[ cl[i] ] = rho_average
            end
            if rho_average > bord_rho[ cl[j] ]
                bord_rho[ cl[j] ] = rho_average
            end
        end
    end
end
for i=1:size_data
    if rho[i] < bord_rho[ cl[i] ]
        halo[i] = 0
    end
end
bord_rho=[]
