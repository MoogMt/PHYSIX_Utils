GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

n_blobs=5
centers=rand(5,2)
nb_points=ones(Int,5)*1000
widths=rand(5)*0.2
points=clustering.createBlobs(nb_points,centers,widths)


file_out=open(string("/home/moogmt/blobs-",n_blobs,"-2.dat"),"w")
for i=1:size(points)[1]
	for j=1:size(points)[2]
		write(file_out,string(points[i,j]," "))
	end
	write(file_out,string("\n"))
end
close(file_out)


# Compute the distance between all points
distance_matrix=clustering.computeDistanceMatrix( points )

#==============================================================================#

# Cluster parameters
n_clusters=5
n_repeat=10
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/cluster-kmenoid-",n_clusters,"-2.dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)


#==============================================================================#

cut_off = 0.02
# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, cluster_indexs = clustering.dauraClustering( distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

file_out=open(string("/home/moogmt/cluster-daura-",cut_off,"-1.dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

#==============================================================================#

n_blobs=5
centers=rand(n_blobs,2)
nb_points=ones(Int,n_blobs)*500
widths=rand(n_blobs)*0.2
points=clustering.createBlobs(nb_points,centers,widths)

file_out=open(string("/home/moogmt/blobs-",n_blobs,"-DPC.dat"),"w")
for i=1:size(points)[1]
	for j=1:size(points)[2]
		write(file_out,string(points[i,j]," "))
	end
	write(file_out,string("\n"))
end
close(file_out)


distance_matrix=clustering.computeDistanceMatrix( points )

cut_off=0.05
file_out=string("/home/moogmt/decision-diagram-",cut_off,".dat")
rho,delta,index,nearest_neighbor=clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix , cut_off, file_out )

min_rho=60
min_delta=0.015
cluster_index,cluster_centers=clustering.densityPeakClusteringSecondStep( rho, delta, index, nearest_neighbor, min_rho, min_delta )

file_out=open(string("/home/moogmt/cluster-DPC-dc-",cut_off,"-",min_rho,"-",min_delta,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_index[i],"\n"))
end
close(file_out)

file_out=open(string("/home/moogmt/cluster-dc-",cut_off,"-centers.dat"),"w")
for i=1:size(cluster_centers)[1]
	write(file_out,string(points[cluster_centers[i],1]," ",points[cluster_centers[i],2],"\n"))
end
close(file_out)
