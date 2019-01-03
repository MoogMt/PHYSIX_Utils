GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"cubefile.jl"))
include(string(CO2folder,"markovCO2.jl"))

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

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
	distance_matrix=computeDistanceMatrix( data , n_dim, max_v, min_v )

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

	# Computing maximum distance between two points
	max_distance=max_array(distance_matrix)

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
	distance_matrix=computeDistanceMatrix( data , n_dim, max_v, min_v )

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

	# Computing maximum distance between two points
	max_distance=max_array(distance_matrix)

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


# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1 # 0.1% cut-off for a state to be considered statististically viable

min_lag=1    # min tau
max_lag=5001 # max tau
d_lag=5      # delta tau
unit=0.005   # units of the simulation

V=8.82
nb_steps=250

distance_data=[]
elf_data=[]
density_data=[]
n_dim=5
distance_configurations=zeros(Int((nbC+nbO)*nb_steps),n_dim)
for step=1:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")
    atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
    density = cube_mod.readCube( string(folder_base,step,"_density.cube") )[3]
    cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
    atoms=cell_mod.wrap(atoms,cell)
    for i=1:nbC
		distances_sort=zeros(Real,nbO+nbC)
		elf_sort=zeros(Real,nbO+nbC)
        for j=i+1:nbC+nbO
			distances_sort[j]=cell_mod.distance(atoms,cell,i,j)
			elf_sort[j] = cube_mod.dataInTheMiddleWME( atoms, cell , i, j, elf ))
            # if cell_mod.distance(atoms,cell,i,j) < 2.5
            #     push!(distance_data,cell_mod.distance(atoms,cell,i,j))
            #     push!(elf_data,cube_mod.dataInTheMiddleWME( atoms, cell , i, j, elf ))
			# 	push!(density_data,cube_mod.dataInTheMiddleWME( atoms, cell , i, j, density ))
            # end
        end
		for j=1:nbO
			for k=j+1:nbO
				if distances_sort[j] > distances_sort[k]
					stock=distances_sort[j]
					distances_sort[j]=distances_sort[k]
					distances_sort[k]=stock
					stock=elf_sort[j]
					elf_sort[j]=elf_sort[k]
					elf_sort[k]=stock
				end
			end
		end
		distance_configuration=distances_sort[1:n_dim]
		elf_configuration=elf_sort[1:n_dim]
    end
end

# nb_box=200
# delta_denself=1/nb_box
# min_dist=100
# max_dist=0
# for i=1:size(distance_data)[1]
#     if distance_data[i] < min_dist
#         global min_dist=distance_data[i]
#     end
#     if distance_data[i] > max_dist
#         global max_dist=distance_data[i]
#     end
# end
# delta_dist=(max_dist-min_dist)/nb_box
#
# hist2D=zeros(nb_box,nb_box)
# count_v=0
# for i=1:size(distance_data)[1]
#     for k=1:nb_box
#         if distance_data[i] > min_dist+(k-1)*delta_dist &&  distance_data[i] < min_dist+k*delta_dist
#             for l=1:nb_box
#                 if elf_data[i] > (l-1)*delta_denself && elf_data[i] < l*delta_denself
#                     hist2D[k,l] += 1
#                     break
#                 end
#             end
#             break
#         end
#     end
#     global count_v += 1
# end
#
# hist2D /= count_v
#
# file_out=open(string(folder_base,"All-elf_hist-",nb_steps,".dat"),"w")
# for i=1:nb_box
# 	for j=1:nb_box
# 		write(file_out,string(min_dist+i*delta_dist," ",j*delta_denself," ",hist2D[i,j],"\n"))
# 	end
# 	write(file_out,string("\n"))
# end
# close(file_out)
#
# hist2D=zeros(nb_box,nb_box)
# count_v=0
# for i=1:size(distance_data)[1]
#     for k=1:nb_box
#         if distance_data[i] > min_dist+(k-1)*delta_dist &&  distance_data[i] < min_dist+k*delta_dist
#             for l=1:nb_box
#                 if density_data[i] > (l-1)*delta_denself && density_data[i] < l*delta_denself
#                     hist2D[k,l] += 1
#                     break
#                 end
#             end
#             break
#         end
#     end
#     global count_v += 1
# end

# hist2D /= count_v
#
# file_out=open(string(folder_base,"All-density_hist-",nb_steps,".dat"),"w")
# for i=1:nb_box
# 	for j=1:nb_box
# 		write(file_out,string(min_dist+i*delta_dist," ",j*delta_denself," ",hist2D[i,j],"\n"))
# 	end
# 	write(file_out,string("\n"))
# end
# close(file_out)
