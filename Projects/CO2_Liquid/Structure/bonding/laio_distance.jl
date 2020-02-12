#==========#
# Laio Algo
#==============================================================================#
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
nb_data_points=Int((nb_steps-1)/1000)
data_set=zeros(Int(nbC*nbO*nb_data_points))
count_=1
for step=1:nb_data_points
    print("Progress:", step/nb_data_points*100,"%\n")
	for carbon=1:nbC
		for oxygen=1:nbO
			data_set[ count_ ] = cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
			global count_ +=1
		end
	end
end

size_data=10000
data_train=data_set[1:size_data,:]

max_v = data_train[1]
min_v = data_train[1]
for i=1:size(data_train)[1]
	if  max_v < data_train[i]
		global max_v = data_train[i]
	end
	if min_v > data_set[i]
		global min_v = data_train[i]
	end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , max_v, min_v )

dc=0.05

rho = gaussianKernel( distance_matrix, dc)

# Sort the rho by decreasing order
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


file_out=open(string(folder,"bond_rho_delta-",dc,"-",size_data,".dat"),"w")
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
    file=open(string(folder,"bond-cluster-",i,"-index-",dc,"-",size_data,".dat"),"w")
    for j=1:size_data
        if i == cl[j]
            write(file,string(data_train[ j ],"\n",))
        end
    end
    close(file)
end
