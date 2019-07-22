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
		data_distance[i,j] = parse(Float64,split(lines2[i])[j+1])
		data_elf[i,j] = parse(Float64,split(lines[i])[j+1])
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
cut_off_dist=0.1
cut_off_all=0.1

file_elf=string(folder_base,"Clust_ELF_decision-diagram-",cut_off_elf,".dat")
file_dist=string(folder_base,"Clust_distance_decision-diagram-",cut_off_dist,".dat")
file_all=string(folder_base,"Clust_all_decision-diagram-",cut_off_all,".dat")

rho_elf,  delta_elf,  index_elf,  nearest_neighbor_elf  = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_elf , cut_off_elf, file_elf )
rho_dist, delta_dist, index_dist, nearest_neighbor_dist = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_distance , cut_off_dist, file_dist )
rho_all,  delta_all,  index_all,  nearest_neighbor_all  = clustering.densityPeakClusteringFirstStepDistanceMatrix( distance_matrix_all , cut_off_all, file_all )


min_rho_elf=10
min_delta_elf=0.4

min_rho_distance=10
min_delta_distance=0.4

cluster_index_elf, cluster_centers_elf = clustering.densityPeakClusteringSecondStep( rho_elf,  delta_elf,  index_elf,  nearest_neighbor_elf,  min_rho_elf,      min_delta_elf )
cluster_index_distance, cluster_centers_distance = clustering.densityPeakClusteringSecondStep( rho_dist, delta_dist, index_dist, nearest_neighbor_dist, min_rho_distance, min_delta_distance )

file_out_elf=open(string(folder_base,"DPC-ELF-dc-",cut_off_elf,".dat"),"w")
file_out_distance=open(string(folder_base,"DPC-distance-dc-",cut_off_distance,".dat"),"w")
file_out_all=open(string(folder_base,"DPC-all-dc-",cut_off_distance,".dat"),"w")
for i=1:nb_point
	for j=1:4
		write(file_out_elf,string(data_elf[i,j]," ") )
		write(file_out_distance,string(data_distance[i,j]," ") )
		write(file_out_all,string(data_distance[i,j]," ") )
	end
	for j=1:4
		write(file_out_all,string(data_elf[i,j]," ") )
	end
	write(file_out_elf,string(cluster_index_elf[i],"\n"))
	write(file_out_distance,string(cluster_index_distance[i],"\n"))
	write(file_out_all,string(cluster_index_elf[i]," ",cluster_index_distance[i],"\n"))
end
close(file_out_elf)
close(file_out_distance)
close(file_out_all)

#==============================================================================#

cut_off_dist=0.1
file_in=open(string(folder_base,"DPC-distance-dc-",cut_off_dist,".dat"))
lines=readlines(file_in)
close(file_in)

nb_lines=size(lines)[1]
nb_dim=size(split(lines[1]))[1]

data=zeros(nb_lines,nb_dim)
for i=1:nb_lines
	for j=1:nb_dim
		data[i,j] = parse(Float64,split(lines[i])[j])
	end
end

cut_off_elf=0.75
cut_off_distance=1.75

nb_err_total=0
nb_err_type1=0
nb_err_type2=0
total_point=0
for i=1:nb_lines
	for j=1:4
		if  data[i,j] > cut_off_elf && data[i,j+4] > cut_off_distance
			nb_err_type1 += 1
			nb_err_total += 1
		elseif data[i,j] < cut_off_elf && data[i,j+4] < cut_off_distance
			nb_err_type2 += 1
			nb_err_total += 1
		end
		total_point += 1
	end
end

file_out_err=open(string(folder_base,"DPC-distance-dc-",cut_off_dist,".dat"),"w")
close(file_out_err)
