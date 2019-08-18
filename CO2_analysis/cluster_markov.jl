GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering


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
