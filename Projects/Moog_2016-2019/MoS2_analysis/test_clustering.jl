
# Laio Algo
GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))

folder2=string("/media/moogmt/Stock/MoS2/MetaExplore/CPMD/Mo3S6/run_2/")

n_dim=9

file_energy_in=open(string(folder2,"ENERGIES"))
lines=readlines(file_energy_in)
close(file_energy_in)
nb_steps=size(lines)[1]
energy=zeros(nb_steps)
energy2=zeros(nb_steps)
for i=1:nb_steps
    energy[i] = parse(Float64,split(lines[i])[4])*(-1)
end

file_colvar_in=open(string(folder2,"COLVAR"))
lines=readlines(file_colvar_in)
close(file_colvar_in)
nb_steps2=size(lines)[1]
colvar=zeros(nb_steps2,n_dim)
count2=0
for i=1:nb_steps2
    if split(lines[i])[1] == "#!"
        continue
    end
    global count2 += 1
    for j=1:n_dim
        colvar[count2,j] = parse(Float64,split(lines[i])[j+1])
    end
end
colvar=colvar[1:2:count2,:]

max_distance=0
for i=1:n_points
    for j=i+1:n_points
        if i==j
            continue
        end
        distance=0
        for k=1:n_dim
            distance += (colvar[i,k]-colvar[j,k])*(colvar[i,k]-colvar[j,k])
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
        if i == j || energy[i] > energy[j]
            continue
        end
        distance=0
        for k=1:n_dim
            distance += (colvar[i,k]-colvar[j,k])*(colvar[i,k]-colvar[j,k])
        end
        if distance < delta[i]
            delta[i] = distance
            nneigh[i] = j
        end
    end
end

max_delta=0
max_energy=0
for i=1:n_points
    if delta[i] > max_delta
        global max_delta = delta[i]
    end
    if energy[i] > max_energy
        global max_energy = energy[i]
    end
end

file_out=open(string(folder2,"decision_diagram.dat"),"w")
for i=1:n_points
    write(file_out,string(energy[i]," ",delta[i],"\n"))
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

n_var=0.1

# min_rho is to be though of carefully, used if there is a group of points
# that are further away than the other, may form an artificial cluster
min_rho=0.5
# Statistically relevant clusters
min_delta=avg_delta+n_var*var_delta

n_cluster=0 # Number of clusters
icl=[]      # Index of cluster centers
cl=ones(Int,n_points)*(-1) # Assignement of the data points (cluster #)

# Determine the cluster centers
for i=1:n_points
    if  delta[i] > 0.5 && max_energy-energy[i] < 0.2
        global n_cluster += 1
        global cl[i] = n_cluster
        global icl=push!(icl,i)
    end
end

print(n_cluster,"\n")

file_out=open(string(folder2,"cluster_center_predicted.dat"),"w")
for i=1:n_cluster
    for j=1:n_dim
        write(file_out,string(colvar[icl[i],j]," "))
    end
    write(file_out,string(energy[icl[i]],"\n"))
end
close(file_out)

file_out_points=open(string(folder2,"points_with_energy.dat"),"w")
for i=1:size(colvar_copy)[1]
    if energy_copy[i] < 170
        continue
    end
    for j=1:n_dim
        write(file_out_points,string(colvar[i,j]," "))
    end
    write(file_out_points,string(energy[i],"\n"))
end
close(file_out_points)
#end
