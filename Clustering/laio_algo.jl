#==============================================================================#
#==========#
# Laio Algo
GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
sim_stride = 1
unit=ps2fs*timestep*sim_stride
nbC=32
nbO=2*nbC

#folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
folder=string("/home/moogmt/CO2/CO2_AIMD/8.82/3000K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.readFastFile( string(folder,file))
cell = cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
n_dim=4
nb_steps_train=400
size_data=nb_steps_train*nbC
data_train=zeros(size_data,n_dim)
for step=1:nb_steps_train
    print("Progress:", step/nb_steps_train*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
		for i=1:n_dim
			data_train[ carbon+nbC*(step-1), i ] = distances[ i ]
		end
	end
end

max_v=data_train[1,:]
min_v=data_train[1,:]
for i=1:size_data
    for j=1:n_dim
        if max_v[j] < data_train[i,j]
            max_v[j] = data_train[i,j]
        end
        if min_v[j] > data_train[i,j]
            min_v[j] = data_train[i,j]
        end
    end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim, max_v, min_v )

dc=0.005

rho = gaussianKernel( distance_matrix, dc)
# Sort the rho by decreasing order
#rho_sorted, index_rho = sort_descend(rho)

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
		global max_rho = rho[i]
	end
end


# put the delta of the max rho point to the max of delta
for i=1:size_data
	if distance_matrix[index[1],i] > delta[1]
		delta[1] = distance_matrix[index[1],i]
	end
end

max_delta=0
for i=1:size(delta)[1]
	if delta[i] > max_delta
		global max_delta  = delta[i]
	end
end

file_out=open(string(folder,"decision-C-",dc,"-",size_data,".dat"),"w")
for i=1:size(rho)[1]
	write(file_out,string(rho[i]/max_rho," ",delta[i]/max_delta,"\n"))
end
close(file_out)

min_delta=0.1
min_rho=0.1

n_cluster=0
cl=ones(Int,size_data)*(-1)
icl=[]
# Determine the cluster centers
for i=1:size_data
    if rho[i] > min_rho && delta[i] > min_delta
		global n_cluster += 1
        cl[index[i]] = n_cluster
        global icl=push!(icl,index[i])
    end
end

for i=1:size_data
    if cl[index[i]] == -1
        cl[index[i]] = cl[ index[nneigh[i]]  ]
    end
end
for i=1:n_cluster
    file=open(string(folder,"coord_predict-C-",dc,"-",size_data,"-cluster-",i,"-index.dat"),"w")
    for j=1:size_data
        if i == cl[j]
            for k=1:n_dim
                write(file,string(data_train[ j, k ]," ",))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end

coord_laio=zeros(n_cluster)
for i=1:n_cluster
	for j=1:n_dim
		if data_train[icl[i],j] < 1.75
			coord_laio[i] += 1
		end
	end
end

agree=0
disagree_lp=0
disagree_lm=0
count_v=1
for step=1:nb_steps_train
	for carbon=1:nbC
		coord_co=0
		for oxygen=1:nbO
			if cell_mod.distance(traj[step],cell,carbon,oxygen+nbC) < 1.75
				coord_co += 1
			end
		end
		if coord_co > coord_laio[cl[count_v]]
			print("STOP - step: ",step," carbon: ",carbon," laio says: ",coord_laio[cl[count_v]]," cut_off says: ",coord_co,"\n")
			global disagree_lm += 1
		elseif coord_co < coord_laio[cl[count_v]]
			print("STOP - step: ",step," carbon: ",carbon," laio says: ",coord_laio[cl[count_v]]," cut_off says: ",coord_co,"\n")
			global disagree_lp += 1
		else
			global agree += 1
		end
		global count_v += 1
	end
end

# TO BE DONE:
# Handshake with ELF (in progress)
# Prediction
# Handshake with O / comparison
# Markovianiation of the results
# Relation between transition and environement
# Function and Fortranization of the code

# Handshake with ELF
include(string(GPfolder,"cubefile.jl"))

folder2=string("/home/moogmt/CO2/CO2_AIMD/ELF/ELF_8.82_results/")

elf_data=[]
elf_distance=[]
n_points=1000
for i=1:n_points
	print("Progress: ",i/n_points*100,"%\n")
	traj, cell_matrix, elf = cube_mod.readCube( string(folder2,i,"_elf.cube"))
	for carbon=1:nbC
		position1=traj.positions[carbon,:]
		for k=1:3
			position1[k]=cell_mod.wrap(position1[k],V)
		end
		for oxygen=1:nbO
			if cell_mod.distance(traj,cell,carbon,nbC+oxygen) < 2.3
				position2=traj.positions[nbC+oxygen,:]
				for j=1:3
				    di = position1[j]-position2[j]
				    if di > V*0.5
						position2[j] += V
				    end
				    if di < -V*0.5
						position2[j] -= V
				    end
				end
				center=(position1+position2)/2.
				for dim=1:3
					center[dim]=cell_mod.wrap(center[dim], V)
				end
				index=[0,0,0]
				for i=1:3
					check=center[i]*elf.nb_vox[i]/V
				    	index[i]=trunc(check)
					if check - index[i] > 0.5
						index[i] += 1
					end
					if index[i] > elf.nb_vox[i]-1
						index[i]=0
					end
				end
				distance1=0
				for k=1:3
				    distance1+=cell_mod.dist1D( index[k]*V/elf.nb_vox[k], center[k], V )^2
				end
				for l=-1:1:1
					for m=-1:1:1
						for n=-1:1:1
							new_index=[0,0,0]
							new_index[1]=index[1]+l
							new_index[2]=index[2]+m
							new_index[3]=index[3]+n
							for l=1:3
						    		if new_index[l] < 0
									new_index[l] = elf.nb_vox[l] - new_index[l]
						    		end
						    		if new_index[l] >= elf.nb_vox[l]-1
									new_index[l] = new_index[l]-elf.nb_vox[l]
						    		end
							end
							distance2=0
							for k=1:3
						    		distance2+=cell_mod.dist1D( new_index[k]*V/elf.nb_vox[k], center[k], V )^2
							end
							if distance2 < distance1
						    		distance1 = distance2
						    		index=new_index
							end
					    	end
					end
				end
				for k=1:3
					index[k] += 1
				end
				push!(elf_data,elf.matrix[index[1],index[2],index[3]])
				push!(elf_distance,cell_mod.distance(traj,cell,carbon,nbC+oxygen))
			end
		end
	end
end

min_=100
max_=0
for i=1:size(elf_distance)[1]
	if min_ > elf_distance[i]
		global min_ = elf_distance[i]
	end
	if max_ < elf_distance[i]
		global max_ = elf_distance[i]
	end
end

n_box=150
delta_distance=(max_-min_)/n_box
delta_elf = 1/n_box

hist2d=zeros(n_box,n_box)
for k=1:size(elf_distance)[1]
	for i=1:n_box
		if elf_distance[k] > min_+(i-1)*delta_distance && elf_distance[k] < min_+i*delta_distance
			for j=1:n_box
				if elf_data[k] > (j-1)*delta_elf && elf_data[k] < j*delta_elf
					hist2d[i,j] += 1
					break
				end
			end
		end
	end
end

file_out=open(string(folder2,"elf_hist.dat"),"w")
for i=1:n_box
	for j=1:n_box
		write(file_out,string(min_+i*delta_distance," ",j*delta_elf," ",hist2d[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)

# Clustering

n_data=10000
data_train=zeros(n_data,2)
data_train[:,1]=elf_distance[1:n_data]
data_train[:,2]=elf_data[1:n_data]

n_dim=2

max_v=data_train[1,:]
min_v=data_train[1,:]
for i=1:n_data
    for j=1:n_dim
        if max_v[j] < data_train[i,j]
            max_v[j] = data_train[i,j]
        end
        if min_v[j] > data_train[i,j]
            min_v[j] = data_train[i,j]
        end
    end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim, max_v, min_v )

dc=0.005

rho = gaussianKernel( distance_matrix, dc)
# Sort the rho by decreasing order
#rho_sorted, index_rho = sort_descend(rho)

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

max_distance=max_array(distance_matrix)

# Compute delta
delta=ones(n_data)*max_distance
delta[ 1 ] = -1
nneigh=zeros(Int,n_data)
for i=2:n_data
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
		global max_rho = rho[i]
	end
end


# put the delta of the max rho point to the max of delta
for i=1:n_data
	if distance_matrix[index[1],i] > delta[1]
		delta[1] = distance_matrix[index[1],i]
	end
end

max_delta=0
for i=1:size(delta)[1]
	if delta[i] > max_delta
		global max_delta  = delta[i]
	end
end

file_out=open(string(folder2,"decision-elf_CO-",dc,"-",n_data,".dat"),"w")
for i=1:size(rho)[1]
	write(file_out,string(rho[i]/max_rho," ",delta[i]/max_delta,"\n"))
end
close(file_out)

min_delta=0.1
min_rho=0.1

n_cluster=0
cl=ones(Int,n_data)*(-1)
icl=[]
# Determine the cluster centers
for i=1:n_data
    if rho[i] > min_rho && delta[i] > min_delta
		global n_cluster += 1
        cl[index[i]] = n_cluster
        global icl=push!(icl,index[i])
    end
end

for i=1:n_data
    if cl[index[i]] == -1
        cl[index[i]] = cl[ index[nneigh[i]]  ]
    end
end
for i=1:n_cluster
    file=open(string(folder2,"Predict-ELF-bond-CO-",dc,"-",n_data,"-cluster-",i,"-index.dat"),"w")
    for j=1:n_data
        if i == cl[j]
            for k=1:n_dim
                write(file,string(data_train[ j, k ]," ",))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end




# HALOOOO
#
# # Do not now exacly what this is...
# halo=copy(cl)
#
# bord_rho=zeros(size_data)
# for i=1:size_data-1
#     for j=i+1:size_data
#         if ( cl[i] == cl[j] ) && ( distance_matrix[i,j] <= dc )
#             rho_average = (rho[i]+rho[j])/2
#             if rho_average > bord_rho[ cl[i] ]
#                 bord_rho[ cl[i] ] = rho_average
#             end
#             if rho_average > bord_rho[ cl[j] ]
#                 bord_rho[ cl[j] ] = rho_average
#             end
#         end
#     end
# end
# for i=1:size_data
#     if rho[i] < bord_rho[ cl[i] ]
#         halo[i] = 0
#     end
# end
# bord_rho=[]
