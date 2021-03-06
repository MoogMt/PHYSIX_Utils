
# Laio Algo
GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))

n_clusters=10
n_dim=2

p_clusters=rand(n_clusters,n_dim)
amplitude=ones(n_clusters)
sigma=ones(n_clusters,n_dim)*0.001

#for n_points in [1000,2000,3000,4000,5000,10000]

n_points=8000

file_out=open(string("/home/moogmt/center_clusters.dat"),"w")
for i=1:n_clusters
    for j=1:n_dim
        write(file_out,string(p_clusters[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

r_points=zeros(Real,n_points,n_dim)
rho=zeros(Real,n_points)
for i=1:n_points
    r_points[i,:]=rand(n_dim)
    for j=1:n_clusters
        exp_int=0
        for k=1:n_dim
            exp_int += (r_points[i,k]-p_clusters[j,k])*(r_points[i,k]-p_clusters[j,k])/(2*sigma[j,k])
        end
        rho[i]+=amplitude[j]*exp(-exp_int)
    end
end

# Sorting
for i=1:n_points
    for j=1:n_points
        if rho[i] > rho[j]
            stock = rho[j]
            rho[j] = rho[i]
            rho[i] = stock
            stock = r_points[j,:]
            r_points[j,:] = r_points[i,:]
            r_points[i,:] = stock
        end
    end
end

file_out=open(string("/home/moogmt/points.dat"),"w")
for i=1:n_points
    for j=1:n_dim
        write(file_out,string(r_points[i,j]," "))
    end
    write(file_out,string(rho[i],"\n"))
end
close(file_out)

max_distance=0
for i=1:n_points
    for j=i+1:n_points
        if i==j
            continue
        end
        distance=0
        for k=1:n_dim
            distance += (r_points[i,k]-r_points[j,k])*(r_points[i,k]-r_points[j,k])
        end
        if distance > max_distance
            global max_distance = distance
        end
    end
end

delta=ones(Real,n_points)*max_distance
nneigh=zeros(Int,n_points)

# Slow
for i=1:n_points
    for j=1:n_points
        if i == j || rho[i] > rho[j]
            continue
        end
        distance=0
        for k=1:n_dim
            distance += (r_points[i,k]-r_points[j,k])*(r_points[i,k]-r_points[j,k])
        end
        if distance < delta[i]
            delta[i] = distance
            nneigh[i] = j
        end
    end
end

function recursiveSearch( used::Vector{T1}, nneigh::Vector{T2}, start::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    for i=start+1:n_points
        if used[i] < 0
            continue
        end
        if nneigh[i] == start
            used[i] = -1
            recursiveSearch( used, nneigh, i)
        end
    end
    return
end

# Destroy non usable points
used=zeros(n_points)

recursiveSearch(used,nneigh,1)

file_out=open(string("/home/moogmt/points2.dat"),"w")
for i=1:n_points
    if used[i] < 0
        continue
    end
    for j=1:n_dim
        write(file_out,string(r_points[i,j]," "))
    end
    write(file_out,string(rho[i],"\n"))
end
close(file_out)

max_delta=0
for i=1:n_points
    if used[i] < 0
        continue
    end
    if delta[i] > max_delta
        global max_delta = delta[i]
    end
end
delta[1]=max_delta

file_out=open(string("/home/moogmt/decision_diagram-",n_points,".dat"),"w")
for i=1:n_points
    if used[i] < 0
        continue
    end
    write(file_out,string(rho[i]/rho[1]," ",delta[i],"\n"))
end
close(file_out)

avg_delta=0
var_delta=0
for i=1:n_points
    global avg_delta += delta[i]
    global var_delta += delta[i]*delta[i]
end
avg_delta /= n_points
var_delta = var_delta/n_points - avg_delta*avg_delta

n_var=0.2

# min_rho is to be though of carefully, used if there is a group of points
# that are further away than the other, may form an artificial cluster
min_rho=0.2
# Statistically relevant clusters
min_delta=avg_delta+n_var*var_delta

n_cluster=0 # Number of clusters
icl=[]      # Index of cluster centers
cl=ones(Int,n_points)*(-1) # Assignement of the data points (cluster #)

# Determine the cluster centers
for i=1:n_points
    if rho[i]/max_rho > min_rho && delta[i]/max_delta > min_delta
        global n_cluster += 1
        global cl[i] = n_cluster
        global icl=push!(icl,i)
    end
end

print(n_points," ",n_cluster,"\n")

file_out=open(string("/home/moogmt/cluster_center_predicted-",n_points,".dat"),"w")
for i=1:n_cluster
    for j=1:n_dim
        write(file_out,string(r_points[icl[i],j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)
#end
