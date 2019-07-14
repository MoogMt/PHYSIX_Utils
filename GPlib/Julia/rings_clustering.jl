GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"clustering.jl"))

# TEST K-medoid
#==============================================================================#
# Definition of the points

n_dim=2
center_ring=zeros(n_dim)
n_points=5000
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

# Cluster parameters
n_clusters=10
n_repeat=50

n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/kmenoid-ring-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)


#=============================================================#
cut_off=0.1

cluster_centers, cluster_sizes, index_cluster = clustering.dauraClustering( size(points)[1] , distance_matrix , cut_off )

file_out=open(string("/home/moogmt/kmenoid-ring-",n_clusters,".dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)
