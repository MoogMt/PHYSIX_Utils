GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

file_in=open(string(folder_base,"1-4_ELF.dat"))
lines=readlines(file_in)
close(file_in)


file_in=open(string(folder_base,"1-4_distance.dat"))
lines2=readlines(file_in)
close(file_in)

data_distance = zeros(size(lines)[1],4)
data_elf = zeros(size(lines)[1],4)

nb_point=size(lines)[1]

for i=1:nb_point
	for j=1:4
		data_distance[i,j] = parse(Float64,split(lines2[i])[j])
		data_elf[i,j] = parse(Float64,split(lines2[i])[j])
	end
end

distance_matrix_elf=zeros(nb_point,nb_point)
distance_matrix_distance=zeros(nb_point,nb_point)
distance_matrix_all=zeros(nb_point,nb_point)
for i=1:nb_point
	for j=i+1:nb_point
		dist_elf=0
		dist_distance=0
		dist_total=0
		for k=1:4
			dist_elf += (data_elf[i,k]-data_elf[j,k])*(data_elf[i,k]-data_elf[j,k])
			dist_distance += (data_distance[i,k]-data_distance[j,k])*(data_distance[i,k]-data_distance[j,k])
		end
		distance_matrix_elf[i,j] = sqrt(dist_elf)
		distance_matrix_elf[j,i] = distance_matrix_elf[i,j]
		distance_matrix_distance[i,j] = sqrt(dist_distance)
		distance_matrix_distance[j,i] = distance_matrix_distance[i,j]
		distance_matrix_all[i,j] = sqrt(dist_elf+dist_distance)
		distance_matrix_all[j,i] = distance_matrix_all[i,j]
	end
end

cut_off_elf=0.05
cut_off_dist=0.05
cut_off_all=0.1

file_elf=string(folder_base,"Clust_ELF_decision-diagram_elf-",cut_off_elf,".dat")
file_dist=string(folder_base,"Clust_ELF_decision-diagram_dist-",cut_off_dist,".dat")
file_all=string(folder_base,"Clust_ELF_decision-diagram_all-",cut_off_all,".dat")

rho_elf,  delta_elf,  index_elf,  nearest_neighbor_elf  = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_elf , cut_off_elf, file_elf )
rho_dist, delta_dist, index_dist, nearest_neighbor_dist = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_distance , cut_off_distance, file_dist )
rho_all,  delta_all,  index_all,  nearest_neighbor_all  = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_all , cut_off_all, file_all )






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
