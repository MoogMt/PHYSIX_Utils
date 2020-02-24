GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering


# TEST K-medoid
#==============================================================================#
# Definition of the points


points=rand(10000,2)
nb_points=size(points)[1]
n_dim=size(points)[2]

file_out=open(string("/home/moogmt/sample_uniform.dat"),"w")
for i=1:nb_points
	for j=1:n_dim
		write(file_out,string(points[i,j]," "))
	end
	write(file_out,string("\n"))
end
close(file_out)


# Compute the distance between all points
distance_matrix=clustering.computeDistanceMatrix( points )

# Cluster parameters

n_clusters=5
n_repeat=10
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/tesselation-kmenoid-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

n_clusters=10
n_repeat=10
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/tesselation-kmenoid-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)


n_clusters=20
n_repeat=10
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/tesselation-kmenoid-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

n_clusters=100
n_repeat=10
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/tesselation-kmenoid-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

#==============================================================================#

cut_off = 0.1
# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, cluster_indexs = clustering.dauraClustering( distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

file_out=open(string("/home/moogmt/tesselation-daura-",cut_off,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)


cut_off = 0.05
# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, cluster_indexs = clustering.dauraClustering( distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

file_out=open(string("/home/moogmt/tesselation-daura-",cut_off,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

cut_off = 0.01
# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, cluster_indexs = clustering.dauraClustering( distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

file_out=open(string("/home/moogmt/tesselation-daura-",cut_off,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

cut_off = 0.005
# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, cluster_indexs = clustering.dauraClustering( distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

file_out=open(string("/home/moogmt/tesselation-daura-",cut_off,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)
