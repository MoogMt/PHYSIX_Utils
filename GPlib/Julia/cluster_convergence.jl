# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=10
co=cut_off[1]

nbC=32
nbO=2*nbC

function voronoiAssign( n_clusters, n_structures, distance_matrix, cluster_sizes, cluster2i, i2cluster )
    cluster_sizes=zeros(n_clusters)
    for i=1:n_structures
        dmin=1.0

    end
    return
end

function swap(table, index1, index2 )
    stock=table[index1]
    table[index1]=table[index2]
    table[index2]=table[stock]
    return
end

function kmenoid_clustering{ T1 <: Int, T2 <: Int, T3 <: Real }( n_clusters::T1 , n_structures::T2, distance_matrix::Vector{Real} )

    cluster_centers=zeros(n_clusters)
    prob=ones(n_structures)
    available=ones(n_structures)

    # Initial Choice

    # initial random choice
    cluster_centers[1]=Int(trunc(rand()*n_structure))+1
    available[cluster_centers[1]] = 0
    for i=2:n_clusters
        # Compute normalized distances to the other cluster centers
        for k=1:n_structures
            dmin=1.0
            for j=1:i-1
                dtmp=distance_matrix[k,cluster_centers[j]]
                if dtmp < dmin
                    dtmp=dmin
                end
            end
            proba[k] = dmin^2
        end
        dtmp=sum(proba)
        proba/=dtmp
        for k=2:nb_steps
            proba[k] += proba[k-1]
        end
        # Look for new cluster center
        found=false
        while ! found
            r=rand()
            if r < proba[1] && available[1] != 0
                cluster_center[i]=1
                good[1]=0
                found = true
            else
                for k=2:n_structures
                    if available[k] != 0 &&  r > proba[k-1] && r < proba[k]
                        cluster_centers[i]=k
                        good[k]=0
                        found=true
                        break
                    end
                end
            end
        end
    end

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
            ass[]
        end
    end

    return cluster_centers
end
