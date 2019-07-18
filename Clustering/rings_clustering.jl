GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

# TEST K-medoid
#==============================================================================#
# Definition of the points

n_dim=2
center_ring=zeros(n_dim)
n_points=10000
radius=2
width=0.2

points=zeros(n_points,n_dim)
points[1:Int(n_points/2),:]=clustering.createRing(Int(n_points/2),center_ring,radius,width)
radius=1
width=0.3
points[Int(n_points/2)+1:n_points,:]=clustering.createRing(Int(n_points/2),center_ring,radius,width)
file_out=open(string("/home/moogmt/ring.dat"),"w")
for i=1:n_points
	for j=1:n_dim
		write(file_out,string(points[i,j]," "))
	end
	write(file_out,string("\n"))
end

# Compute the distance between all points
distance_matrix=clustering.computeDistanceMatrix( points )

# K-Menoid
#=============================================================#
# Cluster parameters
n_clusters=10

for n_clusters in n_clusters_list
	n_repeat=50

	n_structures=size(points)[1]
	cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

	file_out=open(string("/home/moogmt/kmenoid-ring-",n_clusters,".dat"),"w")
	for i=1:size(points)[1]
		write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
	end
	close(file_out)

	file_out=open(string("/home/moogmt/kmenoid-ring-dc-",cut_off,"-centers.dat"),"w")
	for i=1:size(cluster_centers)[1]
		write(file_out,string(points[cluster_centers[i],1]," ",points[cluster_centers[i],2],"\n"))
	end
	close(file_out)

end

# DAURA
#=============================================================#
cut_offs=[0.05,0.1,0.5,1.0,2.0]

for cut_off in cut_offs

	cluster_centers, cluster_sizes, index_cluster = clustering.dauraClustering( distance_matrix , cut_off )

	file_out=open(string("/home/moogmt/daura-ring-dc-",cut_off,".dat"),"w")
	for i=1:size(points)[1]
		write(file_out,string(points[i,1]," ",points[i,2]," ",index_cluster[i],"\n"))
	end
	close(file_out)

	file_out=open(string("/home/moogmt/daura-ring-dc-",cut_off,"-centers.dat"),"w")
	for i=1:size(cluster_centers)[1]
		write(file_out,string(points[cluster_centers[i],1]," ",points[cluster_centers[i],2],"\n"))
	end
	close(file_out)

end
#=============================================================#

cut_off=0.1
file_out=string("/home/moogmt/decision-diagram",cut_off,".dat")

rho,delta,index,nearest_neighbor=clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix , cut_off, file_out )

min_rho=180
min_delta=3
cluster_index,cluster_centers=clustering.densityPeakClusteringSecondStep( rho, delta, index, nearest_neighbor, min_rho, min_delta )

file_out=open(string("/home/moogmt/DPC-ring-dc-",cut_off,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_index[i],"\n"))
end
close(file_out)

file_out=open(string("/home/moogmt/DPC-ring-dc-",cut_off,"-centers.dat"),"w")
for i=1:size(cluster_centers)[1]
	write(file_out,string(points[cluster_centers[i],1]," ",points[cluster_centers[i],2],"\n"))
end
close(file_out)
