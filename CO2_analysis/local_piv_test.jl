GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"clustering.jl"))
include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))


# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1
min_lag=1
max_lag=5001
d_lag=5
unit=0.005


T=3000
V=8.82

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")
#folder_out=string(folder_in)

nbC=32
nbO=64
nb_atoms=nbC+nbO

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)



d=1.75
n=30
m=60

nb_type=2
step_max=100
max_neigh=5
states_sig=zeros(max_neigh*nb_type,nbC,step_max)
#states_sig_sec=zeros(max_neigh*nb_type*max_neigh*nb_type,nbC,step_max)
for step_sim=1:step_max
    print("progress: ",step_sim/step_max*100,"%\n")
    bond_matrix=zeros(Real,96,96)
    index = utils.sequenceMatrixH(nb_atoms)
    for atom1=1:nb_atoms
        for atom2=atom1+1:nb_atoms
            value=utils.switchingFunction(cell_mod.distance(traj[step_sim],cell,atom1,atom2),d,n,m)
            bond_matrix[atom1,atom2]=value
            bond_matrix[atom2,atom1]=value
        end
        # Sorting
        for carbon1=1:nbC
            for carbon2=carbon1+1:nbC
                if bond_matrix[atom1,carbon1] < bond_matrix[atom1,carbon2]
                    stock=bond_matrix[atom1,carbon2]
                    bond_matrix[atom1,carbon2] = bond_matrix[atom1,carbon1]
                    bond_matrix[atom1,carbon1] =stock
                    stock = index[atom1,carbon2]
                    index[atom1,carbon2] = index[atom1,carbon1]
                    index[atom1,carbon1] = stock
                end
            end
        end
        for oxygen1=1:nbO
            for oxygen2=oxygen1+1:nbO
                if bond_matrix[atom1,nbC+oxygen1] < bond_matrix[atom1,nbC+oxygen2]
                    stock=bond_matrix[atom1,nbC+oxygen2]
                    bond_matrix[atom1,nbC+oxygen2] = bond_matrix[atom1,nbC+oxygen1]
                    bond_matrix[atom1,nbC+oxygen1] = stock
                    stock = index[atom1,nbC+oxygen2]
                    index[atom1,nbC+oxygen2] = index[atom1,nbC+oxygen1]
                    index[atom1,nbC+oxygen1] = stock
                end
            end
        end

    end

    for carbon=1:nbC
        states_sig[1:max_neigh,carbon,step_sim]=bond_matrix[carbon,1:max_neigh]
        states_sig[max_neigh+1:2*max_neigh,carbon,step_sim]=bond_matrix[carbon,nbC+1:nbC+max_neigh]
    end
end

nb_points=step_max*nbC

dc=0.05

dist_max=0
file_out=open(string(folder_base,"dist_histogram_Cstates.dat"),"w")
rho=zeros(nb_points)
for step=1:step_max
    print("Progress: ",step/(step_max)*100,"%\n")
    for carbon=1:nbC
        for step2=step+1:step_max
            for carbon2=carbon+1:nbC
                dist=0
                for i=1:size(states_sig)[1]
                    dist += (states_sig[i,carbon,step]-states_sig[i,carbon2,step2])*(states_sig[i,carbon,step]-states_sig[i,carbon2,step2])
                end
                dist=sqrt(dist)
                write(file_out,string(dist,"\n"))
                if dist_max < dist
                    global dist_max = dist
                end
                rho[(step-1)*nbC+carbon] += exp(-dist/dc)
                rho[(step2-1)*nbC+carbon2] += exp(-dist/dc)
            end
        end
    end
end
close(file_out)
#
# file_in=open(string(folder_base,"dist_histogram_Cstates.dat"))
# lines=readlines(file_in)
# close(file_in)
# nb_points2=size(lines)[1]
# dist=zeros(nb_points2)
# for i=1:nb_points2
#     dist[i] = parse(Float64,lines[i])
# end
# min_dist=0
# max_dist=dist_max
# nb_box=100
# delta=(max_dist-min_dist)/nb_box
# hist1D=zeros(nb_box)
# for data=1:nb_points2
#     print("Progress: ",data/nb_points2*100,"%\n")
#     for box=1:nb_box
#         if dist[data] > min_dist+(box-1)*delta && dist[data] < min_dist+box*delta
#             hist1D[box] += 1
#             break
#         end
#     end
# end
# hist1D /= sum(hist1D)
# file_out=open(string(folder_base,"hist_dist.dat"),"w")
# for i=1:nb_box
#     write(file_out,string(min_dist+i*delta," ",hist1D[i],"\n"))
# end
# close(file_out)

delta=ones(nb_points)*dist_max
nneigh=zeros(Int,step_max,nbC)
for step=1:step_max
    print("Progress: ",step/(step_max)*100,"%\n")
    for carbon=1:nbC
        for step2=1:step_max
            for carbon2=1:nbC
                if carbon==carbon2 && step==step2
                    continue
                end
                if rho[(step-1)*nbC+carbon] > rho[(step2-1)*nbC+carbon2]
                    continue
                end
                dist=0
                for i=1:size(states_sig)[1]
                    dist += (states_sig[i,carbon,step]-states_sig[i,carbon2,step2])*(states_sig[i,carbon,step]-states_sig[i,carbon2,step2])
                end
                dist=sqrt(dist)
                if dist < delta[(step-1)*nbC+carbon]
                    delta[(step-1)*nbC+carbon] = dist
                    nneigh[step,carbon]=(step2-1)*nbC+carbon2
                end
            end
        end
    end
end

file_out=open(string(folder_base,"decision_diagram_markov.dat"),"w")
for i=1:nb_points
    write(file_out,string(rho[i]," ",delta[i],"\n"))
end
close(file_out)


min_delta=0.18
n_cluster=0 # Number of clusters
icl=[]      # Index of cluster centers
cl=ones(Int,nb_points)*(-1) # Assignement of the data points (cluster #)

# Determine the cluster centers
for i=1:nb_points
    if  delta[i] > min_delta
        global n_cluster += 1
        global cl[i] = n_cluster
        global icl=push!(icl,i)
    end
end

for i=1:nb_points
    if cl[i] == -1
        cl[i] = cl[ nneigh[i] ]
    end
end


for i=1:n_cluster
    file_out=open(string(folder_base,"cluster-",i,".dat"),"w")
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
