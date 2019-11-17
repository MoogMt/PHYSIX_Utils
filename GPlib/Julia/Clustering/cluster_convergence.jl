GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering


# TEST K-medoid
#==============================================================================#
# Definition of the points

centers=ones(2,2)
centers[1,:]=zeros(2)

points=clustering.createBlobs([200,200],centers,[0.02,0.1])
nb_points=size(points)[1]
n_dim=size(points)[2]

file_out=open(string("/home/moogmt/blobs.dat"),"w")
for i=1:nb_points
	for j=1:n_dim
		write(file_out,string(points[i,j]," "))
	end
	write(file_out,string("\n"))
end
close(file_out)


# Compute the distance between all points
distance_matrix=clustering.computeDistanceMatrix( points )

# Cluster parameters
n_clusters=2

n_repeat=100
n_structures=size(points)[1]
cluster_indexs, cluster_centers, cluster_sizes = clustering.kmedoidClustering( n_structures, distance_matrix, n_clusters, n_repeat )

file_out=open(string("/home/moogmt/test-kmenoid-1.dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",cluster_indexs[i],"\n"))
end
close(file_out)

points2=rand(2000,2)+ones(2000,2)*0.4
assignments2=clustering.voronoiAssign( points, n_clusters, cluster_centers, points2 )

file_out=open(string("/home/moogmt/test-kmenoid-2.dat"),"w")
for i=1:size(points)[1]
	write(file_out,string(points[i,1]," ",points[i,2]," ",assignments2[i],"\n"))
end
close(file_out)


# TEST Daura
#==============================================================================#
# Definition of the points
ndim=2
offset=1.2
n_points=2000
points=zeros(n_points,ndim)
for i=1:n_points
    points[i,:]=rand(ndim)
end
points[Int(trunc(n_points/2))+1:n_points,:] += offset

max=points[1,:]
min=points[1,:]
for i=1:size(points)[1]
    for j=1:2
        if max[j] < points[i,j]
            max[j] = points[i,j]
        end
        if min[j] > points[i,j]
            min[j] = points[i,j]
        end
    end
end

# Compute the distance between all points
distance_matrix=computeDistanceMatrix( points)

cut_off = 0.8

# Cluster parameters
n_structures=size(points)[1]
cluster_centers, cluster_sizes, index_cluster = dauraClustering( size(points)[1] , distance_matrix , cut_off )
n_clusters=size(cluster_centers)[1]

n_test=200
points2=vcat(rand(n_test,2), rand(n_test,2)+offset )
points2_assignment=voronoiAssign( points, n_clusters, cluster_centers, points2 , max, min )

using PyPlot

figure()
for i=1:size(points2)[1]
	if points2_assignment[i] == 1
		plot( points2[i,1] , points2[i,2] , "r."  )
	elseif points2_assignment[i] == 2
		plot( points2[i,1] , points2[i,2] , "b."  )
    elseif points2_assignment[i] == 3
        plot( points2[i,1] , points2[i,2] , "g."  )
    elseif points2_assignment[i] == 4
        plot( points2[i,1] , points2[i,2] , "c."  )
	end
end
for i=1:n_clusters
    plot( points[cluster_centers[ i ] , 1] , points[cluster_centers[ i ], 2] , "kd"  )
end

#==========#
# K-medoid
#==============================================================================#

include("contactmatrix.jl")

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
stride = 1
unit=ps2fs*timestep*stride
start_time=5
start_step=Int(start_time/unit)
nbC=32
nbO=2*nbC

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.read( string(folder,file), stride, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
maxN=4
nb_dim=maxN+1
data_set=zeros(nb_steps*nbC,nb_dim)
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

n_train=4000
n_dim_analysis=9
data_train=data_set[1:n_train,1:n_dim_analysis ]
data_predict=data_set[n_train+1:nb_steps,1:n_dim_analysis ]

max=data_train[1,:]
min=data_train[1,:]
for i=1:n_train
    for j=1:n_dim_analysis
        if max[j] < data_train[i,j]
            max[j] = data_train[i,j]
        end
        if min[j] > data_train[i,j]
            min[j] = data_train[i,j]
        end
    end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim_analysis, max, min )

# Cluster parameters
n_clusters=3
precision=0.00000000000001
n_repeat=10

print("Clustering\n")
cluster_indexs, cluster_centers, cluster_sizes, assignments = kmedoidClustering( n_train, distance_matrix, n_clusters, precision , n_repeat)

print("Printing cluster centers\n")
file = open( string( folder, "center_cluster-",n_dim_analysis,".dat" ), "w" )
for i=1:n_clusters
    for j=1:n_dim_analysis
        write( file, string(data_train[ cluster_centers[i] , j ]," " ) )
    end
    write( file, "\n" )
end
close( file )

print("Printing Training Clusters\n")
for i=1:size(assignments)[1]
    file=open( string( folder, string("cluster",i,"-",n_dim_analysis,".dat") ), "w" )
    for j=1:size(assignments)[2]
        for k=1:n_dim_analysis
            if assignments[ i, j ] != 0
                write(file,string( data_train[ assignments[ i, j ], k ] ," ") )
            end
        end
        write(file,string("\n"))
    end
    close(file)
end

# Voronoi assignment
print("Predicting assignments\n")
predict_assignment = voronoiAssign( data_train, n_clusters, cluster_centers,  data_predict, max, min )

print("Printing Assignments\n")
for cluster=1:n_clusters
    file=open( string( folder, string("cluster",cluster,"-",n_dim_analysis,".dat") ), false,true,false,false,true )
    for elements=1:nb_steps-n_train
        if predict_assignment[ elements] == cluster
            for data=1:n_dim_analysis
                write(file,string( data_predict[ elements , data ]," "))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end

#Oxygen
#===#

include("contactmatrix.jl")

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
stride = 1
unit=ps2fs*timestep*stride
start_time=5
start_step=Int(start_time/unit)
nbC=32
nbO=2*nbC

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.read( string(folder,file), stride, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
nb_dim=3
data_set=zeros(nb_steps*nbO,nb_dim)
fileO=open(string(folder,"distancesNN_O.dat"),"w")
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for oxygen=1:nbO
		distances = zeros(nbC)
		for carbon=1:nbC
			distances[carbon] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
        maxN=2
		for i=1:maxN
			data_set[ oxygen+nbO*(step-1), i ] = distances[ i ]
            write(fileO, string(distances[i]," ") )
		end
        count=maxN+1
        for i=1:maxN-1
            for j=i+1:maxN
                a=cell_mod.distance(traj[step],cell,Int(index[i]),oxygen)
                b=cell_mod.distance(traj[step],cell,Int(index[j]),oxygen)
                c=cell_mod.distance(traj[step],cell,Int(index[i]),Int(index[j]))
                angle=acosd((a*a+b*b-c*c)/(2*a*b))
                write(fileO,string(angle," "))
                data_set[ oxygen+nbO*(step-1), count ] = angle
                count += 1
            end
        end
        write(fileO,string("\n"))
	end
end
close(fileO)

n_train=100
n_dim_analysis=3
data_train=data_set[1:n_train,1:n_dim_analysis ]
data_predict=data_set[n_train+1:nb_steps,1:n_dim_analysis ]
data_set=[]

max=data_train[1,:]
min=data_train[1,:]
for i=1:n_train
    for j=1:n_dim_analysis
        if max[j] < data_train[i,j]
            max[j] = data_train[i,j]
        end
        if min[j] > data_train[i,j]
            min[j] = data_train[i,j]
        end
    end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim_analysis, max, min )

# Cluster parameters
n_clusters=2
precision=0.00000000000001
n_repeat=1

print("Clustering\n")
cluster_indexs, cluster_centers, cluster_sizes, assignments = kmedoidClustering( n_train, distance_matrix, n_clusters, precision , n_repeat)

print("Printing cluster centers\n")
file = open( string( folder, "center_clusterO-",n_dim_analysis,".dat" ), "w" )
for i=1:n_clusters
    for j=1:n_dim_analysis
        write( file, string(data_train[ cluster_centers[i] , j ]," " ) )
    end
    write( file, "\n" )
end
close( file )

print("Printing Training Clusters\n")
for i=1:size(assignments)[1]
    file=open( string( folder, string("clusterO",i,"-",n_dim_analysis,".dat") ), "w" )
    for j=1:size(assignments)[2]
        for k=1:n_dim_analysis
            if assignments[ i, j ] != 0
                write(file,string( data_train[ assignments[ i, j ], k ] ," ") )
            end
        end
        write(file,string("\n"))
    end
    close(file)
end

# Voronoi assignment
print("Predicting assignments\n")
predict_assignment = voronoiAssign( data_train, n_clusters, cluster_centers,  data_predict, max, min )

print("Printing Assignments\n")
for cluster=1:n_clusters
    file=open( string( folder, string("clusterO-",cluster,"-",n_dim_analysis,".dat") ), false,true,false,false,true )
    for elements=1:nb_steps-n_train
        if predict_assignment[ elements] == cluster
            for data=1:n_dim_analysis
                write(file,string( data_predict[ elements , data ]," "))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end

#========#
# Daura
#==============================================================================#

include("contactmatrix.jl")

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
stride = 1
unit=ps2fs*timestep*stride
start_time=5
start_step=Int(start_time/unit)
nbC=32
nbO=2*nbC

folder=string("/home/moogmt/CO2/CO2_AIMD/",V,"/",T,"K/")
file="TRAJEC_wrapped.xyz"

traj = filexyz.read( string(folder,file), stride, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

# Training set
nb_dim=10
data_set=zeros(nb_steps*nbC,nb_dim)
fileC=open(string(folder,"distancesNN.dat"),"w")
for step=1:nb_steps
    print("Progress:", step/nb_steps*100,"%\n")
	for carbon=1:nbC
		distances = zeros(nbO)
		for oxygen=1:nbO
			distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
		end
		index=sortmatrix( distances )
        maxN=4
		for i=1:maxN
			data_set[ carbon+nbC*(step-1), i ] = distances[ i ]
            write(fileC, string(distances[i]," ") )
		end
        a=cell_mod.distance(traj[step],cell,carbon,Int(index[1]+nbC))
        b=cell_mod.distance(traj[step],cell,carbon,Int(index[2]+nbC))
        c=cell_mod.distance(traj[step],cell,Int(index[1]+nbC),Int(index[2]+nbC))
        angle=acosd((a*a+b*b-c*c)/(2*a*b))
        count=maxN+1
        for i=1:maxN-1
            for j=i+1:maxN
                a=cell_mod.distance(traj[step],cell,carbon,Int(index[i]+nbC))
                b=cell_mod.distance(traj[step],cell,carbon,Int(index[j]+nbC))
                c=cell_mod.distance(traj[step],cell,Int(index[i]+nbC),Int(index[j]+nbC))
                angle=acosd((a*a+b*b-c*c)/(2*a*b))
                write(fileC,string(angle," "))
                data_set[ carbon+nbC*(step-1), count ] = angle
                count += 1
            end
        end
        write(fileC,string("\n"))
	end
end
close(fileC)

n_train=4000
n_dim_analysis=5
data_train=data_set[1:n_train,1:n_dim_analysis ]
data_predict=data_set[n_train+1:nb_steps,1:n_dim_analysis ]
data_set=[]

max=data_train[1,:]
min=data_train[1,:]
for i=1:n_train
    for j=1:n_dim_analysis
        if max[j] < data_train[i,j]
            max[j] = data_train[i,j]
        end
        if min[j] > data_train[i,j]
            min[j] = data_train[i,j]
        end
    end
end

print("Computing distance matrix\n")
distance_matrix=computeDistanceMatrix( data_train , n_dim_analysis, max, min )

cut_off = 0.6

print("Clustering\n")
cluster_centers, cluster_sizes, index_cluster = dauraClustering( n_train , distance_matrix , cut_off )

n_clusters=size(cluster_centers)[1]

print("Printing cluster centers\n")
file = open( string( folder, "center_cluster-",n_dim_analysis,"_data_",cut_off,".dat" ), "w" )
for i=1:n_clusters
    for j=1:n_dim_analysis
        write( file, string(data_train[ cluster_centers[i] , j ]," " ) )
    end
    write( file, "\n" )
end
close( file )


n=3
file=open(string(folder,"cluster-",n,"-",n_dim_analysis,"_daura_",cut_off,".dat"),"w")
for j=1:n_train
    if n == index_cluster[j]
        for k=1:n_dim_analysis
            write(file, string( data_train[ j  , k ] ," ") )
        end
        write(file,string("\n"))
    end
end
close(file)
