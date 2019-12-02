GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Counts the percent of C with at least three O neighbors
# For all Temperatures and Volumes
# Needs wrapped trajectory files with thermalization
# steps removed

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz

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
size_data=6000
n_dim=4
data_train=zeros(nb_steps*nbC,n_dim)
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		for i=1:maxN
			data_train[ carbon+nbC*(step-1), i ] = distances[ i ]
		end
		index=sortmatrix( distances )
		for i=1:maxN
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

file_out=open(string(folder,"test_test.dat"),"w")
for i=1:size(rho)[1]
	write(file_out,string(rho[i]/max_rho," ",delta[i],"\n"))
end
close(file_out)

min_delta=0.1
min_rho=0.1*max_rho

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
    file=open(string(folder,"cluster-",i,"-index.dat"),"w")
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


# Do not know exacly what this is...
halo=copy(cl)

bord_rho=zeros(size_data)
for i=1:size_data-1
    for j=i+1:size_data
        if ( cl[i] == cl[j] ) && ( distance_matrix[i,j] <= dc )
            rho_average = (rho[i]+rho[j])/2
            if rho_average > bord_rho[ cl[i] ]
                bord_rho[ cl[i] ] = rho_average
            end
            if rho_average > bord_rho[ cl[j] ]
                bord_rho[ cl[j] ] = rho_average
            end
        end
    end
end
for i=1:size_data
    if rho[i] < bord_rho[ cl[i] ]
        halo[i] = 0
    end
end
bord_rho=[]
