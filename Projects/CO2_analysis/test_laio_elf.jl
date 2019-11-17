
# Laio Algo
GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))
include(string(GPfolder,"cubefile.jl"))

folder2=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/")
traj, cell_matrix, elf = cube_mod.readCube( string(folder2,100,"_elf.cube"))

nb_vox=size(elf.matrix)[1]

max_distance=nb_vox*sqrt(3)

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



max_rho=0
max_delta=0
for i=1:n_points
    if rho[i] > max_rho
        max_rho = rho[i]
    end
    if delta[i] > max_delta
        max_delta = delta[i]
    end
end

file_out=open(string("/home/moogmt/decision_diagram-",n_points,".dat"),"w")
for i=1:n_points
    write(file_out,string(rho[i]/max_rho," ",delta[i],"\n"))
end
close(file_out)

avg_delta=0
var_delta=0
for i=1:n_points
    avg_delta += delta[i]
    var_delta += delta[i]*delta[i]
end
avg_delta /= n_points
var_delta = var_delta/n_points - avg_delta*avg_delta

n_var=3

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
        n_cluster += 1
        cl[i] = n_cluster
        icl=push!(icl,i)
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
