# TESTED AND FUNCTIONNALS
#==============================================================================#
function computeDistance{ T1 <: Real, T2 <: Real, T3 <: Int }( data::Array{T1}, data_point::Vector{T2}, index::T3 )
	return sum( ( data[ index, :  ] - data_point[:] ).*( data[ index, :  ] - data_point[:] ) )
end
function computeDistance{ T1 <: Real, T2 <: Int, T3 <: Int }( data::Array{T1}, index1::T2, index2::T3 )
	return sum( ( data[ index1, : ]-data[ index2 , : ] ).*( data[ index1, : ]-data[ index2 , : ] ) )
end
function computeDistance{ T1 <: Real,T2 <: Int,  T3 <: Real, T4 <: Real, T5 <: Int, T6 <: Int }( data::Array{T1}, n_dim::T2 , max::Vector{T3}, min::Vector{T4}, index1::T5, index2::T6 )
    dist=0
    for i=1:n_dim
        dist += ( ((data[index1,i]-min[i])/(max[i]-min[i])) - ((data[index2,i]-min[i])/(max[i]-min[i])) )^2
    end
	return dist
end
function computeDistanceMatrix{ T1 <: Real }( data::Array{T1} )
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, i, j )
        end
    end
    return matrix
end
function computeDistanceMatrix{ T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }( data::Array{T1}, n_dim::T2, max::Vector{T3}, min::Vector{T4} )
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, n_dim, max, min, i, j )
        end
    end
    return matrix
end
function swap{ T1 <: Real, T2 <: Int, T3 <: Int }( table::Vector{T1}, index1::T2, index2::T3 )
    stock=table[index1]
    table[index1]=table[index2]
    table[index2]=stock
    return
end
function initializeCenters{ T1 <: Int, T2 <: Real , T3 <: Int}( n_structures::T1 , distance_matrix::Array{T2,2}  , n_clusters::T3  )
    # Bookkeep
    available=ones(Int,n_structures)

    cluster_centers=zeros(Int,n_clusters)
    prob=zeros(n_structures)

    cluster_centers[1]= trunc( rand()*n_structures )   + 1
    available[ cluster_centers[1] ] = 0
    for i=2:n_clusters
        for k=1:n_structures
            min_distance=1
            for j=1:i-1
                dist_temp=distance_matrix[ k, Int(cluster_centers[j]) ]
                if dist_temp < min_distance
                    min_distance = dist_temp
                end
            end
            prob[k] = min_distance^2
        end
        prob=cumsum(prob/sum(prob))
        found=false
        while ! found
            r=rand()
            if r < prob[1] && available[1] == 1
                cluster_centers[i] = 1
                available[1] = 0
                found = true
            else
                for k=2:n_structures
                    if r > prob[k-1] && r < prob[k] && available[k] == 1
                        cluster_centers[i] = k
                        available[k] = 0
                        found = true
                        break
                    end
                end
            end
        end
    end

    return cluster_centers
end
function voronoiAssign{ T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real }( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_point::Vector{T4} )
	index_cluster=1
	min_dist=sum( ( data[ cluster_centers[1], :  ] - data_point ).*( data[ cluster_centers[1], :  ] - data_point ) )
	for i=2:n_clusters
		dist=sum( ( data[ cluster_centers[i], :  ] - data_point ).*( data[ cluster_centers[i], :  ] - data_point ) )
		if dist < min_dist
			min_dist=dist
			index_cluster=i
		end
	end
	return index_cluster
end
function voronoiAssign{ T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real }( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_points::Array{T4} )
	n_points=size(data_points)[1]
	point_clusters=zeros(Int,n_points)
	for i=1:n_points
		point_clusters[i] = voronoiAssign( data, n_clusters, cluster_centers, data_points[i,:] )
	end
	return point_clusters
end
function voronoiAssignSingle{ T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int}( distance_matrix::Array{T1,2}, nb_clusters::T2, cluster_centers::Vector{T3}, index::T4 )
    min_dist=distance_matrix[ index ,cluster_centers[1] ]
    index_cluster=1
    for i=2:nb_clusters
    dist=distance_matrix[ index ,cluster_centers[i] ]
        if dist < min_dist
            min_dist = dist
            index_cluster = i
        end
    end
    return index_cluster
end
function voronoiAssignAll{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int}( nb_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4} )
    # Index of the cluster for each structure
    cluster_indexs = zeros(Int, nb_structures )
    assignments = zeros(Int, nb_clusters, nb_structures )
    # Contains sizer of clusters
    cluster_sizes=zeros(Int,nb_clusters)
    for structure=1:nb_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[structure], cluster_sizes[ cluster_indexs[structure] ] ] = structure
    end
    return cluster_indexs, cluster_sizes, assignments
end
function voronoiAssignAll{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Int, T7 <: Int }( n_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4}, cluster_indexs::Vector{T5}, cluster_sizes::Vector{T6}, assignments::Array{T7,2} )
    cluster_sizes=zeros(Int,nb_clusters)
    for structure=1:n_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[ structure ], cluster_sizes[ cluster_indexs[ structure ] ] ] = structure
    end
    return
end
function computeCost{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int }( n_structures::T1, distance_matrix::Array{T2,2}, cluster_centers::Vector{T3} , cluster_indexs::Vector{T4} )
    cost=0
    for i=1:n_structures
        cost += distance_matrix[ i, cluster_centers[ cluster_indexs[i] ] ]
    end
    return cost
end
function updateCenters{ T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int , T5 <: Int }( distance_matrix::Array{T1,2} , n_clusters::T2, cluster_centers::Vector{T3}, cluster_sizes::Vector{T4}, assignments::Array{T5,2} )
    for cluster=1:n_clusters
        new_center = cluster_centers[ cluster ]
        cost_min = sum(distance_matrix[ assignments[ cluster,1], assignments[ cluster ,1:cluster_sizes[ cluster ] ] ])
        for structure_in_cluster=2:cluster_sizes[ cluster ]
            cost = sum(distance_matrix[ assignments[ cluster, structure_in_cluster ], assignments[ cluster ,1:cluster_sizes[ cluster ] ] ])
            if cost < cost_min
                cost_min=cost
                new_center = assignments[ cluster , structure_in_cluster ]
            end
        end
        cluster_centers[ cluster ] = new_center
    end
end
function sortmatrix( x )
    sizex=size(x)[1]
    index_x=zeros(sizex)
    for i=1:sizex
        index_x[i]=i
    end
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stocki=index_x[i]
                x[i]=x[j]
                index_x[i]=index_x[j]
                x[j]=stock
                index_x[j]=stocki
            end
        end
    end
    return index_x
end
function kmedoidClustering{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Real }( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3 , precision::T4 )
    # Initialization of centers
    cluster_centers=initializeCenters(n_structures, distance_matrix, n_clusters )
    # Assign all clusters
    cluster_indexs, cluster_sizes, assignments =voronoiAssignAll( n_structures, distance_matrix, n_clusters, cluster_centers )
    # Compute original cost
    old_cost=computeCost(n_structures,distance_matrix,cluster_centers,cluster_indexs)

    old_cost=1
    while true
        voronoiAssignAll( n_structures, distance_matrix, n_clusters, cluster_centers,cluster_indexs, cluster_sizes, assignments  )
        cost=computeCost( n_structures, distance_matrix, cluster_centers, cluster_indexs)
        if abs(cost-old_cost) < precision
            break
        end
        old_cost=cost
        updateCenters( distance_matrix, n_clusters, cluster_centers, cluster_sizes, assignments )
    end
    return cluster_indexs, cluster_centers, cluster_sizes, assignments
end
function kmedoidClustering{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Int }( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3 , precision::T4, n_repeat::T5 )
    # First Kmenoid
    cluster_indexs_best, cluster_centers_best, cluster_sizes_best, assignments_best = kmedoidClustering( n_structures, distance_matrix, n_clusters, precision )
    old_cost=computeCost( n_structures, distance_matrix, cluster_centers_best, cluster_indexs_best )
    for i=2:n_repeat
        cluster_indexs, cluster_centers, cluster_sizes, assignments = kmedoidClustering( n_structures, distance_matrix, n_clusters, precision )
        cost=computeCost( n_structures, distance_matrix, cluster_centers, cluster_indexs )
        if old_cost > cost
            cluster_indexs_best = cluster_indexs
            cluster_centers_best = cluster_centers
            cluster_sizes_best = cluster_sizes
            assignments_best = assignments
        end
    end
    return cluster_indexs_best, cluster_centers_best, cluster_sizes_best, assignments_best
end
#==============================================================================#


# TEST
#==============================================================================#
# Definition of the points
points=zeros(200,2)
for i=1:50
    points[i,:]=rand(2)
end
for i=51:100
    points[i,:]=rand(2)
end
points[51:100,1] += 0.8
for i=101:150
    points[i,:]=rand(2)
end
points[101:150,2] += 0.8
for i=151:200
    points[i,:]=rand(2)
end
points[151:200,1] += 0.8
points[151:200,2] += 0.8

# Compute the distance between all points
distance_matrix=computeDistanceMatrix( points )

# Cluster parameters
n_clusters=4
precision=0.00000000000001
n_structures=size(points)[1]

cluster_indexs, cluster_centers, cluster_sizes, assignments = kmedoidClustering( n_structures, distance_matrix, n_clusters, precision )
old_cost=computeCost( n_structures, distance_matrix, cluster_centers, cluster_indexs )

points2=rand(2000,2)+0.4
points2_assignment=voronoiAssign( points, n_clusters, cluster_centers, points2 )

using PyPlot

figure(1)
plot( points[cluster_centers[ 1 ] ,1 ] , points[cluster_centers[ 1 ],2] , "rd"  )
plot( points[cluster_centers[ 2 ] ,1 ] , points[cluster_centers[ 2 ],2] , "bd"  )
plot( points[cluster_centers[ 3 ] ,1 ] , points[cluster_centers[ 3 ],2] , "gd"  )
plot( points[cluster_centers[ 4 ] ,1 ] , points[cluster_centers[ 4 ],2] , "kd"  )
for i=1:size(points2)[1]
	if points2_assignment[i] == 1
		plot( points2[i,1] , points2[i,2] , "r."  )
	elseif points2_assignment[i] == 2
		plot( points2[i,1] , points2[i,2] , "b."  )
	elseif points2_assignment[i] == 3
		plot( points2[i,1] , points2[i,2] , "g."  )
	elseif points2_assignment[i] == 4
		plot( points2[i,1] , points2[i,2] , "k."  )
	end
end
for j=1:cluster_sizes[ 1 ]
    plot( points[ assignments[ 1, j ] , 1 ] , points[ assignments[ 1, j ] ,2 ] , "r."  )
end
for j=1:cluster_sizes[ 2 ]
    plot( points[ assignments[ 2, j ] , 1 ] , points[ assignments[ 2, j ] ,2 ] , "b."  )
end
for j=1:cluster_sizes[ 3 ]
    plot( points[ assignments[ 3, j ] , 1 ] , points[ assignments[ 3, j ] ,2 ] , "g."  )
end
for j=1:cluster_sizes[ 4 ]
    plot( points[ assignments[ 4, j ] , 1 ] , points[ assignments[ 4, j ] ,2 ] , "k."  )
end

#==============================================================================#

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

folder=string("/home/moogmt/CO2/CO2_AIMD/",V,"/",T,"K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.read( string(folder,file), stride, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
nb_dim=10
data_set=zeros(nb_steps*nbC,nb_dim)
max=zeros(nb_dim)
min=ones(nb_dim)*100000
fileC=open(string(folder,"distancesNN.dat"),"w")
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
        maxN=4
		for i=1:maxN
			data_set[ carbon+nbC*(step-1), i ] = distances[ i ]
            if data_set[ carbon+nbC*(step-1), i ] > max[ i ]
                max[ i ] = data_set[carbon+nbC*(step-1), i ]
            end
            if data_set[ carbon+nbC*(step-1), i ] < min[ i ]
                min[ i ] = data_set[ carbon+nbC*(step-1), i ]
            end
            write(fileC, string(distances[i]," ") )
		end
        a=cell_mod.distance(traj[step],cell,carbon,Int(index[1]+nbC))
        b=cell_mod.distance(traj[step],cell,carbon,Int(index[2]+nbC))
        c=cell_mod.distance(traj[step],cell,Int(index[1]+nbC),Int(index[2]+nbC))
        angle=acosd((a*a+b*b-c*c)/(2*a*b))
        count=maxN+1
        for i=1:maxN-1
            for j=i+1:maxN
                a=cell_mod.distance(traj[step],cell,carbon,Int(index[i]+nbC))
                b=cell_mod.distance(traj[step],cell,carbon,Int(index[j]+nbC))
                c=cell_mod.distance(traj[step],cell,Int(index[i]+nbC),Int(index[j]+nbC))
                angle=acosd((a*a+b*b-c*c)/(2*a*b))
                write(fileC,string(angle," "))
                data_set[ carbon+nbC*(step-1), count ] = angle
                if data_set[ carbon+nbC*(step-1), count ] > max[ count ]
                    max[ count ] = data_set[ carbon+nbC*(step-1), count ]
                end
                if data_set[ carbon+nbC*(step-1), count ] < min[ count ]
                    min[ count ] = data_set[ carbon+nbC*(step-1), count ]
                end
                count += 1
            end
        end
        write(fileC,string("\n"))
	end
end
close(fileC)

n_train=4000
n_dim_analysis=4
data_train=data_set[1:n_train,1:n_dim_analysis]
data_predict=data_set[n_train+1:nb_steps,1:n_dim_analysis]
data_set=[]

distance_matrix=computeDistanceMatrix( data_train , n_dim_analysis, max[1:n_dim_analysis], min[1:n_dim_analysis] )

# Cluster parameters
n_clusters=3
precision=0.00000000000001
n_repeat=10

cluster_indexs, cluster_centers, cluster_sizes, assignments = kmedoidClustering( n_train, distance_matrix, n_clusters, precision , n_repeat)

file = open( string( folder, "center_cluster.dat" ), "w" )
for i=1:n_clusters
    for j=1:n_dim_analysis
        write( file, string(data_train[ cluster_centers[i] , j ]," " ) )
    end
    write( file, "\n" )
end
close( file )

for i=1:size(assignments)[1]
    file=open( string( folder, string("cluster",i,".dat") ), "w" )
    for j=1:size(assignments)[2]
        for k=1:n_dim_analysis
            if assignments[ i, j ] != 0
                write(file,string( data_train[ assignments[ i, j ], k ] ," ") )
            end
        end
        write(file,string("\n"))
    end
    close(file)
end

# Voronoi assignment
predict_assignment = voronoiAssign( data_train, n_clusters, cluster_centers,  data_predict )

for cluster=1:n_clusters
    file=open( string( folder, string("cluster",cluster,".dat") ), false,true,false,false,true )
    for elements=1:nb_steps-n_train
        if predict_assignment[ elements] == cluster
            for data=1:n_dim_analysis
                write(file,string( data_predict[ elements , data ]," "))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end
