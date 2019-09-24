GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion

# Folder for data


folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=1001
d_lag=2
unit=0.005

# for V in Volume
#     for T in Temperatures


T=3000
V=8.82

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

# if ! isfile(file)
#     continue
# end

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]
n_dim=4
nbC=32

data=zeros(nbC*nb_steps,n_dim)

for step=1:nb_steps
    for carbon=1:nbC
        distancesCO=zeros(nbO)
        for oxygen=1:nbO
            distancesCO[oxygen]=cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
        end
        sort!(distancesCO)
        data[(step-1)*nbC+carbon,:] = distancesCO[1:n_dim]
    end
end

nb_data=size(data)[1]

nb_box=100
min_=ones(n_dim)*V
max_=zeros(n_dim)
for i=1:nb_data
    for j=1:n_dim
        if min_[j] > data[i,j]
            min_[j] = data[i,j]
        end
        if max_[j] < data[i,j]
            max_[j] = data[i,j]
        end
    end
end

delta=zeros(n_dim)
for i=1:n_dim
    delta[i] = (max_[i]-min_[i])/(nb_box-1)
end

hist_data=zeros(nb_box,n_dim)
for i=1:nb_data
    for j=1:n_dim
        hist_data[ Int(round( (data[i,j]-min_[j])/delta[j] ))+1,j ] += 1
    end
end
for i=1:4
    hist_data[:,i] /= sum(hist_data[:,i])
end

file_out=open(string(folder_out,"distances_data.dat"),"w")
for i=1:nb_box
    for j=1:n_dim
        write(file_out,string(i*delta[j]+min_[j]," ",hist_data[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

nb_box_markov=[10,10]
delta_markov=zeros(2)

delta_markov[1] = (max_[3]-min_[3])/(nb_box_markov[1]-1)
delta_markov[2] = (max_[4]-min_[4])/(nb_box_markov[2]-1)
n_dim_markov=size(nb_box_markov)[1]

states_markov=zeros(nb_data,2)
hist_markov=zeros(nb_box_markov[1],nb_box_markov[2])
sum_=0
for i=1:nb_data
    hist_markov[ Int(round(( data[i,3]-min_[3])/delta_markov[1] ) )+1 , Int(round((data[i,4]-min_[4])/delta_markov[2]))+1 ] +=1
end
hist_markov_percent = hist_markov/sum(hist_markov)

file_out=open(string(folder_out,"map_3-4C-",nb_box_markov[1],"-",nb_box_markov[2],".dat"),"w")
for i=1:nb_box_markov[1]
    for j=1:nb_box_markov[2]
        write(file_out,string(i*delta_markov[1]+min_[3]," ",j*delta_markov[2]+min_[4]," ",hist_markov_percent[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

nb_states=prod(nb_box_markov)
transition_matrix=zeros(nb_states,nb_states)


#==============================================================================#

data=zeros(nb_steps,nbC,n_dim)

for step=1:nb_steps
    for carbon=1:nbC
        distancesCO=zeros(nbO)
        for oxygen=1:nbO
            distancesCO[oxygen]=cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
        end
        sort!(distancesCO)
        data[step,carbon,:] = distancesCO[1:n_dim]
    end
end

nb_box=100
min_=ones(nbC,n_dim)*V
max_=zeros(nbC,n_dim)
for carbon=1:nbC
    for step=1:nb_steps
        for j=1:n_dim
            if min_[carbon,j] > data[step,carbon,j]
                min_[carbon,j] = data[step,carbon,j]
            end
            if max_[carbon,j] < data[step,carbon,j]
                max_[carbon,j] = data[step,carbon,j]
            end
        end
    end
end

delta=zeros(nbC,n_dim)
for carbon=1:nbC
    for i=1:n_dim
        delta[carbon,i] = (max_[carbon,i]-min_[carbon,i])/(nb_box-1)
    end
end

hist_data=zeros(nbC,nb_box,n_dim)
for i=1:nb_data
    for j=1:n_dim
        hist_data[ Int(round( (data[i,j]-min_[j])/delta[j] ))+1,j ] += 1
    end
end
for i=1:4
    hist_data[:,i] /= sum(hist_data[:,i])
end

file_out=open(string(folder_out,"distances_data.dat"),"w")
for i=1:nb_box
    for j=1:n_dim
        write(file_out,string(i*delta[j]+min_[j]," ",hist_data[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

nb_box_markov=[10,10]
delta_markov=zeros(2)

delta_markov[1] = (max_[3]-min_[3])/(nb_box_markov[1]-1)
delta_markov[2] = (max_[4]-min_[4])/(nb_box_markov[2]-1)
n_dim_markov=size(nb_box_markov)[1]

states_markov=zeros(nb_data,2)
hist_markov=zeros(nb_box_markov[1],nb_box_markov[2])
sum_=0
for i=1:nb_data
    hist_markov[ Int(round(( data[i,3]-min_[3])/delta_markov[1] ) )+1 , Int(round((data[i,4]-min_[4])/delta_markov[2]))+1 ] +=1
end
hist_markov_percent = hist_markov/sum(hist_markov)

file_out=open(string(folder_out,"map_3-4C-",nb_box_markov[1],"-",nb_box_markov[2],".dat"),"w")
for i=1:nb_box_markov[1]
    for j=1:nb_box_markov[2]
        write(file_out,string(i*delta_markov[1]+min_[3]," ",j*delta_markov[2]+min_[4]," ",hist_markov_percent[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

nb_states=prod(nb_box_markov)
transition_matrix=zeros(nb_states,nb_states)
