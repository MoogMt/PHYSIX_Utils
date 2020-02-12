file=open("/home/moogmt/test.dat")
lines=readlines(file)
close(file)

nb_line=size(lines)[1]

nb_points = 0
for i=1:nb_line
    line = split(lines[i])
    if isempty(line)
        continue
    else
        global nb_points += 1
    end
end

rho=zeros(Real,nb_points)
delta=zeros(Real,nb_points)
x=zeros(Real,nb_points)
y=zeros(Real,nb_points)

count_point=1
for i=1:nb_line
    line = split(lines[i])
    if isempty(line)
        continue
    else
        x[count_point] = parse(Float64,line[1])
        y[count_point] = parse(Float64,line[2])
        rho[count_point] = parse(Float64,line[3])*(-1)
        global count_point += 1
    end
end

n_clusters_base=3
x_cluster=rand(n_clusters_base,1)
y_cluster=rand(n_clusters_base,1)

nb_points=10000
rho=zeros(Real,nb_points)
delta=zeros(Real,nb_points)
x=zeros(Real,nb_points)
y=zeros(Real,nb_points)
for i=1:nb_points
    x[i]=rand()
    y[i]=rand()
    for j=1:n_clusters_base
        rho[i] += exp( -((x[i]-x_cluster[j])^2+(y[i]-y_cluster[j])^2)/0.005 )
    end
end

file_test=open(string("/home/moogmt/test2.dat"),"w")
for i=1:nb_points
    write(file_test,string(x[i]," ",y[i]," ",rho[i],"\n"))
end
close(file_test)

GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"clustering.jl"))

index=clustering.simpleSequence(nb_points)

for i=1:nb_points
    for j=i+1:nb_points
        if rho[i] < rho[j]
            stock = rho[j]
            rho[i] = rho[j]
            rho[j] = stock
            stock = index[i]
            index[i] = index[j]
            index[j] = stock
        end
    end
end

max_distance=0
for i=1:nb_points
    for j=i+1:nb_points
        dist=sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) )
        if dist > max_distance
            global max_distance = dist
        end
    end
end

delta=ones(nb_points)*max_distance
delta[1] = -1
nneigh=zeros(Int,nb_points)

for i=2:nb_points
    for j=1:i-1
        distance_ij  = ( x[ index[i]] - x[ index[j] ] )*( x[ index[i] ]- x[ index[j] ] )
        distance_ij += ( y[ index[i]] - y[ index[j] ] )*( y[ index[i] ]- y[ index[j] ] )
        distance_ij = sqrt( distance_ij)
        if distance_ij < delta[ i ]
            delta[ i ] = distance_ij
            nneigh[ i ] = j
        end
    end
end

max_rho=rho[1]

max_delta=0
for i=1:nb_points
    if delta[i] > max_delta
        global max_delta = delta[i]
    end
end
delta[1]=max_delta

file_out=open(string("/home/moogmt/CO2.dat"),"w")
for i=1:nb_points
    write(file_out,string(rho[i]/max_rho," ",delta[i]/max_delta,"\n"))
end
close(file_out)

n_cluster=0 # Number of clusters
icl=[]      # Index of cluster centers
cl=ones(Int,nb_points)*(-1) # Assignement of the data points (cluster #)

# Determine the cluster centers
for i=1:nb_points
    if rho[i]/max_rho > 0.5 && delta[i]/max_delta > 0.6
        global n_cluster += 1
        global cl[index[i]] = n_cluster
        global icl=push!(icl,index[i])
    end
end

for i=1:nb_points
    if cl[index[i]] == -1
        cl[index[i]] = cl[ index[nneigh[i]]  ]
    end
end

for i=1:n_cluster
    file_out=open(string("/home/moogmt/cluster-",i,".dat"),"w")
    for j=1:nb_points
        if cl[j] == i
            write(file_out,string(x[j]," ",y[j],"\n"))
        end
    end
    close(file_out)
end

file_out=open(string("/home/moogmt/center_cluster.dat"),"w")
for i=1:n_cluster
    write(file_out,string(x[icl[i]]," ",y[icl[i]],"\n"))
end
close(file_out)
