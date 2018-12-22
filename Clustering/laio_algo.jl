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
maxN=4
n_dim=maxN#+1
data_set=zeros(nb_steps*nbC,n_dim)
fileC=open(string(folder,"distancesNN.dat"),"w")
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
		for i=1:maxN
			data_set[ carbon+nbC*(step-1), i ] = distances[ i ]
            write(fileC, string(distances[i]," ") )
		end
        # a=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[1]))
        # b=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[2]))
        # c=cell_mod.distance(traj[step],cell,nbC+Int(index[1]),nbC+Int(index[2]))
        # angle=acosd((a*a+b*b-c*c)/(2*a*b))
        #data_set[ carbon+nbC*(step-1), maxN+1 ] = angle
        #write(fileC, string(angle," ") )
        write(fileC,string("\n"))
	end
end
close(fileC)

size_data=10000
data_train=data_set[1:size_data,:]

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

dc=0.001
rho=gaussianKernel(distance_matrix,dc)
max_distance=max_array(distance_matrix)
rho_sorted, index_rho = sort_descend(rho)

delta=zeros(size_data)
delta[ index_rho[1] ] = -1
nneigh=zeros(Int,size_data)
for i=2:size_data
    delta[ index_rho[i] ] = max_distance
    for j=1:i-1
        if distance_matrix[ index_rho[i] , index_rho[j] ] < delta[ index_rho[i] ]
            delta[ index_rho[i] ] = distance_matrix[ index_rho[i], index_rho[j] ]
            nneigh[ index_rho[i] ] = index_rho[j]
        end
    end
end

max_rho=0
for i=1:size(rho)[1]
 	if rho[i] > max_rho
		global max_rho = rho[i]
	end

end
delta[ index_rho[1] ] = max_vector( delta )
file_out=open(string(folder,"rho_delta-",dc,"-",size_data,".dat"),"w")
for i=1:size(rho)[1]
	write(file_out,string(rho[i]/max_rho," ",delta[i],"\n"))
end
close(file_out)


min_delta=0.10
min_rho=0.05*max_rho

n_cluster=0
cl=ones(Int,size_data)
cl*=-1
icl=[]
for i=1:size_data
    if rho[i] > min_rho && delta[i] > min_delta
        global n_cluster+=1
        cl[i] = n_cluster
        global icl=push!(icl,i)
    end
end
for i=1:size_data
    if cl[ index_rho[i] ] == -1
        cl[ index_rho[i] ] = cl[ nneigh[ index_rho[i] ] ]
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
