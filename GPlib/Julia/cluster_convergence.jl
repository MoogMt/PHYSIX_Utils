    # Loading file
include("contactmatrix.jl")

# TESTED AND FUNCTIONNALS
function compute_distance( data )
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = sqrt(sum( (points[i,:]-points[j,:]).*(points[i,:]-points[j,:]) ))
        end
    end
    return matrix
end
function swap(table, index1, index2 )
    stock=table[index1]
    table[index1]=table[index2]
    table[index2]=stock
    return
end
function initialize_medoid( distance_matrix , n_structures , n_clusters  )
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
    cluster_sizes = zeros(Int, nb_clusters)
    for structure=1:nb_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[structure], cluster_sizes[ cluster_indexs[structure] ] ] = cluster_indexs[ structure ]
    end
    return cluster_indexs, cluster_sizes, assignments
end
function voronoiAssignAll{ T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Int, T7 <: Int }( nb_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4}, cluster_indexs::Vector{T5}, cluster_sizes::Vector{T6}, assignments::Array{T7,2} )
    for structure=1:nb_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignemnts[ cluster_indexs[structure], cluster_sizes[ cluster_indexs[structure] ] ] = cluster_indexs[ structure ]
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
# Definition of the points
points=zeros(100,3)
for i=1:50
    points[i,:]=rand(3)
end
for i=51:100
    points[i,:]=rand(3)+10
end

# Compute the distance between all points
distance_matrix=compute_distance(points)

# Cluster parameters
nb_clusters=2
nb_structures=size(points)[1]

# Initialization of centers
cluster_centers=initialize_medoid(distance_matrix, nb_structures, nb_clusters )
# Assign all clusters
cluster_indexs, cluster_sizes, assignments =voronoiAssignAll( nb_structures, distance_matrix, nb_clusters, cluster_centers )
# Compute original cost
cost=computeCost(nb_structures,distance_matrix,cluster_centers,cluster_indexs)


function kmenoid_clustering{ T1 <: Int, T2 <: Int, T3 <: Real }( n_clusters::T1 , n_structures::T2, distance_matrix::Vector{Real} )

    cluster_centers=zeros(n_clusters)

    initialize_medoid( distance_matrix, n_structures, n_ clusters, cluster_centers

    old_cost=1

    # Minimization of medoids
    #  - Look up table
    cluster2i=zeros(n_clusters,n_structures)
    i2cluster=zeros(n_structures)
    # - Size of clusters
    cluster_sizes=zeros(n_clusters)
    while true
        # Voronoi assignment
        cluster_sizes=zeros(n_clusters)
        for i=1:n_structures
            dmin=1.0
            for j=1:n_structures
                if distance_matrix[i,j] < dmin
                    dmin = distance_matrix[i,j]
                    i2cluster[i]=j
                end
            end
            j=i2cluster[i]
            cluster_size[j] += 1
            cluster2i[j,cluster_size[j]]=i
        end

        # Cost of the clustering
        cost=0
        for i=1:n_structures
            cost += distance_matrix[i,cluster2i[i]]
        end
        # Exit if nothing changed (convergence reached)
        if cost == old_cost
            break
        end

        # Updating Medoid
        for i=1:n_clusters
            dmin=1.
            for j=1:cluster_size[i]
                dtmp = sum(distance_matrix[ cluster2i[i,j], cluster2i[ i, 1:cluster_size[i] ] ])
                if dtmp < dmin
                    dmin=dtmp
                    newcenter=cluster2i[i,j]
                end
            end
            cluster_center[i] = newcenter
        end
    end

    # Sorting cluster
    cluster_member=cluster2i
    cluster_ranking=zeros(n_clusters)
    for i=1:n_clusters
        cluster_ranking[i]=i
    end
    for i=1:n_clusters-1
        for j=i+1:n_clusters
            if cluster_sizes[i] > cluster_sizes[j]
                swap(cluster_size,i,j)
                swap(cluster_ranking,i,j)
            end
        end
    end

    assigned=zeros(n_structures)
    for i=1:n_clusters
        for j=1:cluster_size[i]
            assigned[cluster_member[i,j]] = assigned[cluster_member[i,j]]+1
        end
    end

    for i=1:n_structure
        if assigned[i] != 1
            print("Oups")
        end
    end

    file=open(string(pwd(),"/log-clusters.dat"),"w")
    for i=1:n_clusters
        write(file,string(i," ",cluster[ i,cluster_sizes(i) ] ) )
    end
    close(file)

    return cluster_centers
end
