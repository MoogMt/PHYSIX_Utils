# TESTED AND FUNCTIONNALS
#==============================================================================#
module clustering

export computeDistance, computeDistanceMatrix, computeDistanceMatrixAndMax, swap
export initializeCenter, voronoiAssign, voronoiAssignSingle, voronoiAssignAll
export computeCost, updateCenters, sortMatrix
export kmedoidClustering, computeClusteringCoefficients
export dauraClustering
export gaussianKernel, maxArray, simpleSequence, maxVector
export densityPeakClusteringTrain

function computeDistance( data::Array{T1}, data_point::Vector{T2}, index::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Int }
	return sum( ( data[ index, :  ] - data_point[:] ).*( data[ index, :  ] - data_point[:] ) )
end
function computeDistance( data::Array{T1}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
	return sum( ( data[ index1, : ]-data[ index2 , : ] ).*( data[ index1, : ]-data[ index2 , : ] ) )
end
function computeDistance( data::Array{T1}, n_dim::T2 , max::Vector{T3}, min::Vector{T4}, index1::T5, index2::T6 ) where { T1 <: Real,T2 <: Int,  T3 <: Real, T4 <: Real, T5 <: Int, T6 <: Int }
    dist=0
    for i=1:n_dim
        dist += ( ((data[index1,i]-min[i])/(max[i]-min[i])) - ((data[index2,i]-min[i])/(max[i]-min[i])) )*( ((data[index1,i]-min[i])/(max[i]-min[i])) - ((data[index2,i]-min[i])/(max[i]-min[i])) )
    end
	return dist
end
function computeDistanceMatrix( data::Array{T1} ) where { T1 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, i, j )
        end
    end
    return matrix
end
function computeDistanceMatrix( data::Array{T1}, max::T2, min::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
			matrix[i,j] = sqrt( ( ((data[i]-min)/(max-min)) - ((data[j]-min)/(max-min)) )*( ((data[i]-min)/(max-min)) - ((data[j]-min)/(max-min)) ) )
        end
    end
    return matrix
end
function computeDistanceMatrix( data::Array{T1}, n_dim::T2, max::Vector{T3}, min::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, n_dim, max, min, i, j )
        end
    end
    return matrix
end
function computeDistanceMatrixAndMax( data::Array{T1}, n_dim::T2, max::Vector{T3}, min::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    size_data=size(data)[1]
	max_matrix=0
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, n_dim, max, min, i, j )
			if matrix[i,j] > max_matrix
				max_matrix = matrix[i,j]
			end
        end
    end
    return matrix, max_matrix
end
function swap( table::Vector{T1}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    stock=table[index1]
    table[index1]=table[index2]
    table[index2]=stock
    return
end
function initializeCenters( n_structures::T1 , distance_matrix::Array{T2,2}  , n_clusters::T3  ) where { T1 <: Int, T2 <: Real , T3 <: Int}
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
function voronoiAssign( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_point::Vector{T4} , max::Vector{T5}, min::Vector{T6} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real , T5 <: Real, T6 <: Real }
	index_cluster=1
	min_dist=sum( ( (data[ cluster_centers[1], :  ]-min)./(max-min) - (data_point-min)./(max-min) ).*( (data[ cluster_centers[1], :  ]-min)./(max-min) - (data_point-min)./(max-min) ) )
	for i=2:n_clusters
		dist=sum( ( (data[ cluster_centers[i], :  ]-min)./(max-min) - (data_point-min)./(max-min) ).*( (data[ cluster_centers[i], :  ]-min)./(max-min) - (data_point-min)./(max-min) ) )
		if dist < min_dist
			min_dist=dist
			index_cluster=i
		end
	end
	return index_cluster
end
function voronoiAssign( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_points::Array{T4} , max::Vector{T5}, min::Vector{T6} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real , T5 <: Real, T6 <: Real }
	n_points=size(data_points)[1]
	point_clusters=zeros(Int,n_points)
	for i=1:n_points
		point_clusters[i] = voronoiAssign( data, n_clusters, cluster_centers, data_points[i,:], max, min )
	end
	return point_clusters
end
function voronoiAssignSingle( distance_matrix::Array{T1,2}, nb_clusters::T2, cluster_centers::Vector{T3}, index::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int}
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
function voronoiAssignAll( nb_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int}
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
function voronoiAssignAll( n_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4}, cluster_indexs::Vector{T5}, cluster_sizes::Vector{T6}, assignments::Array{T7,2} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Int, T7 <: Int }
    cluster_sizes=zeros(Int,nb_clusters)
    for structure=1:n_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[ structure ], cluster_sizes[ cluster_indexs[ structure ] ] ] = structure
    end
    return
end
function computeCost( n_structures::T1, distance_matrix::Array{T2,2}, cluster_centers::Vector{T3} , cluster_indexs::Vector{T4} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int }
    cost=0
    for i=1:n_structures
        cost += distance_matrix[ i, cluster_centers[ cluster_indexs[i] ] ]
    end
    return cost
end
function updateCenters( distance_matrix::Array{T1,2} , n_clusters::T2, cluster_centers::Vector{T3}, cluster_sizes::Vector{T4}, assignments::Array{T5,2} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int , T5 <: Int }
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
function kmedoidClustering( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3 , precision::T4 ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Real }
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
function kmedoidClustering( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3 , precision::T4, n_repeat::T5 ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Int }
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
function computeClusteringCoefficients( distance_matrix::Array{T1,2}, n_clusters::T2 , cluster_sizes::Vector{T3} , assignments::Array{T4,2} , cut_off::T5 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int , T5 <: Real }
    clustering_coefficients=zeros(n_clusters)
    if cut_off <= 0.0
        return
    end
    for cluster=1:n_clusters
        if cluster_sizes[ cluster ] == 1
            clustering_coefficients[ cluster ] = 1.0
        else
            for i=2:cluster_sizes[ cluster ] - 1
                for j=i+1:cluster_sizes[cluster]
                    dist=distance_matrix[ assignments[ cluster, i], assignments[cluster,j] ]
                    if  dist > 0 && dist < cut_off
                        clustering_coefficients += 1
                    end
                end
            end
        end
        clustering_coefficients[ cluster ] /= cluster_sizes[i]*(cluster_sizes[i]-1)/2
    end
    return clustering_coefficients
end
function dauraClustering( n_elements::T1, distance_matrix::Array{T2,2} , cut_off::T3 ) where { T1 <: Int, T2 <: Real, T3 <: Real }
    cluster_sizes=zeros(Int,n_elements)
    cluster_centers=zeros(Int,n_elements)
    used=zeros(Int,n_elements)
    n_element_left = n_elements
    n_clusters = 0
    while sum(used) < n_elements
        n_neighbours = 0
        index_center = 0
        # Choosing most connected unused point as cluster center
        for i=1:n_elements
            if used[i] != 1
                n_local_neighbours = 0
                for j=i+1:n_elements
                    if used[j] != 1
                        if distance_matrix[i,j] < cut_off
                            n_local_neighbours += 1
                        end
                    end
                end
                if n_local_neighbours > n_neighbours
                    n_neighbours = n_local_neighbours
                    index_center = i
                end
            end
        end
        if n_neighbours != 0
            # We have a new cluster
            n_clusters += 1
            n_element_left -= 1
            used[ index_center ] = 1
            cluster_centers[ n_clusters ] = index_center
            cluster_sizes[ n_clusters ] = n_neighbours+1 # +1 because the center counts...
            # Removing the neighborhood of the chosen cluster
            for n=1:n_elements
                if used[n] != 1 && distance_matrix[ index_center, n ] < cut_off
                    n_element_left -= 1
                    used[ n ] = 1
                end
            end
        else
            # For elements that are aloooone
            n_clusters += 1
            n_element_left -= 1
            used[ n_clusters ] = index_center
            cluster_sizes[ n_clusters ] = 1
        end
    end
    # Removing useless memory slots
    cluster_centers = cluster_centers[1:n_clusters]
    cluster_sizes = cluster_sizes[1:n_clusters]
    used=[]
    # Voronoi assignements
    index_data = zeros(n_elements)
    for element=1:n_elements
        min_dist=100000000
        for i=1:n_clusters
            if cluster_centers[i] != element
                dist = distance_matrix[ cluster_centers[i], element ]
                if dist < min_dist
                    index_data[ element ] = i
                    min_dist=dist
                end
            end
        end
    end
    # Voronoi assignement of the points
    return cluster_centers, cluster_sizes, index_data
end
function gaussianKernel( matrix_distance::Array{T1,2} , cut_off_distance::T2) where { T1 <: Real, T2 <: Real }
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
function maxArray( distance_matrix::Array{T1,2} ) where { T1 <: Real }
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
function maxVector( vector::Vector{T1} ) where { T1 <: Real }
    max=0
    for i=1:size(vector)[1]
        if max < vector[i]
            max=vector[i]
        end
    end
    return max
end
function simpleSequence( size_::T1 ) where { T1 <: Int }
    vector=zeros(Int, size_)
    for i=1:size_
        vector[i] = i
    end
    return vector
end
function sortMatrix( x::Vector{T1} ) where { T1 <: Real }
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
function sortDescend( vector::Vector{T1} ) where { T1 <: Real }
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
function sortDescend( data::Array{T1,2}, index::T2 )  where { T1 <: Real , T2 <: Int }
	size_data=size(data)[1]
	n_dim=size(data)[2]
	for i=1:size_data-1
		for j=i+1:size_data
			if data[i,index] < data[j,index]
				stock = data[i,:]
				data[i,:] = data[j,:]
				data[j,:] = stock
			end
		end
	end
	return
end
function densityPeakClusteringTrain( data::Array{T1,2}, dc::T2 ) where { T1 <: Real, T2 <: Real }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers
	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i] > min_rho && delta[i] > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	print("n_cluster: ",n_cluster,"\n")

	if n_cluster != 0
	# Affectation of points to clusters using their nearest neighbor
		for i=1:size_data
	    	if cl[index[i]] == -1
	        	cl[index[i]] = cl[ index[nneigh[i]]  ]
	    	end
		end
		return cl, icl
	else
		return [], []
	end
	
end
function densityPeakClusteringTrain( data::Array{T1,2}, dc::T2 , min_rho::T3, min_delta::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers
	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i] > min_rho && delta[i] > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	# Affectation of points to clusters using their nearest neighbor
	for i=1:size_data
	    if cl[index[i]] == -1
	        cl[index[i]] = cl[ index[nneigh[i]]  ]
	    end
	end

	return cl, icl
end
function densityPeakClusteringFirstStep( data::Array{T1,2}, dc::T2 , file::T3 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	# Writting decision diagram
	file_out=open(file,"w")
	for i=1:size(rho)[1]
		write(file_out,string(rho[i]," ",delta[i],"\n"))
	end
	close(file_out)

	return rho, delta
end
function densityPeakClusteringFirstStep( data::Array{T1,2}, dc::T2 , min_rho::T3, min_delta::T4, file::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: AbstractString }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	# Writting decision diagram
	file_out=open(file,"w")
	for i=1:size(rho)[1]
		write(file_out,string(rho[i]," ",delta[i],"\n"))
	end
	close(file_out)

	return rho, delta
end
function densityPeakClusteringSecondStep( rho::Vector{T1}, delta::Vector{T2}, min_rho::T3, min_delta::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }

	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers
	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i] > min_rho && delta[i] > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	# Affectation of points to clusters using their nearest neighbor
	for i=1:size_data
	    if cl[index[i]] == -1
	        cl[index[i]] = cl[ index[nneigh[i]]  ]
	    end
	end

	return cl, icl
end

end
