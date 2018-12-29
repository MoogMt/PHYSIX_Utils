GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))

#==========#
# Laio Algo
#==============================================================================#

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
unit=ps2fs*timestep
nbC=32
nbO=2*nbC

#folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
folder=string("/home/moogmt/CO2/CO2_AIMD/",V,"/",T,"K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.readFastFile( string(folder,file))
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
maxN=4
n_dim=maxN
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
        a=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[1]))
        b=cell_mod.distance(traj[step],cell,carbon,nbC+Int(index[2]))
        c=cell_mod.distance(traj[step],cell,nbC+Int(index[1]),nbC+Int(index[2]))
        angle=acosd((a*a+b*b-c*c)/(2*a*b))
        data_set[ carbon+nbC*(step-1), maxN+1 ] = angle
        write(fileC, string(angle," ") )
        write(fileC,string("\n"))
	end
end
close(fileC)

size_data=4000
data_train=data_set[1:size_data,:]

max_v = data_train[1,:]
min_v = data_train[1,:]

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

# n_dim=2
# data_train=zeros(size_data,n_dim)
# points=rand(2000,n_dim)
# points2=rand(2000,n_dim)
# points2=points2+0.5
# data_train[1:2000,:]=points[:,:]
# data_train[2001:4000,:]=points2
#
# PyPlot.figure()
# PyPlot.plot(data_train[:,1],data_train[:,2],"r.")
# PyPlot.show()
#
# max=data_train[1,:]
# min=data_train[1,:]

# for i=1:size_data
#     for j=1:n_dim
#         if max[j] < data_train[i,j]
#             max[j] = data_train[i,j]
#         end
#         if min[j] > data_train[i,j]
#             min[j] = data_train[i,j]
#         end
#     end
# end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim, max_v, min_v )


dc=0.001

rho=zeros(size_data)
for i=1:size_data
	for j=1:size_data
		if distance_matrix[i,j] < dc && i != j
			rho[i] += 1
		end
	end
end

delta=zeros(size_data)
index_delta=zeros(size_data)
for i=1:size_data
	for j=1:size_data
		if i != j && rho[j] > rho[i]
			if delta[i] == 0
				delta[i] = distance_matrix[i,j]
				index_delta[i] = j
			elseif delta[i] > distance_matrix[i,j]
				delta[i] = distance_matrix[i,j]
				index_delta[i] = j
			end
		end
	end
end

index_max=1
max_rho=rho[1]
for i=2:size_data
	if rho[i] > max_rho
		global max_rho = rho[i]
		global index_max=i
	end
end

delta[index_max]=distance_matrix[index_max,1]
for i=2:size_data
	if delta[index_max] < distance_matrix[index_max,i]
		delta[index_max] = distance_matrix[index_max,i]
	end
end



# file=open(string(folder,"data_cluster_test.dat"),"w")
# for i=1:size_data
#     write(file,string(rho[i]," ",delta[i],"\n"))
# end
# close(file)
