include("contactmatrix.jl")

temperature=3000
volume=8.82

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",volume,"/",temperature,"K/")
file_in=string(folder,"TRAJEC_wrapped.xyz")

stride_sim=5
unit=0.0005# in ps

stride_analysis=1
start_time=5
start_step=Int(start_time/(unit*stride_sim))

atoms = filexyz.read( file_in, stride_analysis, start_step )
cell=cell_mod.Cell_param( volume, volume, volume )

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

function searchGroupMember{ T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }( matrix::Array{T1}, list::Vector{T2}, index::T3 , group_nb::T4 )
    for i=1:size(matrix)[1]
        if matrix[index,i] > 0
            if list[i] == 0
                list[i]=group_nb
                list=searchGroupMember(matrix,list,i,group_nb)
            end
        end
    end
    return list
end

cut_off = 1.7

sizes=[]
times=[]
for step=1:nb_steps
    print("Building molecules: ", step/nb_steps*100,"%\n")

    #=================#
    # Building matrix
    #========================================================#
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        for j=i+1:nb_atoms
            if cell_mod.distance(atoms[step],cell,i,j) < cut_off
                matrix[i,j]=1
                matrix[j,i]=1
            end
        end
    end
    #=========================================================#

    #====================================#
    # Exploring the matrix for molecules
    #=========================================================#
    nb_mol=0
    mol_index=zeros(nb_atoms)
    for i=1:nb_atoms
        if mol_index[i] == 0
            nb_mol += 1
            mol_index[i]=nb_mol
            mol_index = searchGroupMember(matrix,mol_index,i,nb_mol)
        end
    end
    #=========================================================#

    #======================#
    # Processing molecules
    #=========================================================#
    for i=1:nb_mol
        size=0
        for j=1:nb_atoms
            if mol_index[j] == i
                size += 1
            end
        end
        push!(times,step)
        push!(sizes,size)
    end
    #=========================================================#

end

sim_time_max=nb_steps*unit*stride_sim

#Making 5ps windows size
time_window=1 # in ps
nb_window=Int(trunc(sim_time_max/time_window)+1)
times *= unit*stride_sim

hist2d=zeros(nb_window,nb_atoms)
nb_sizes=size(sizes)[1]
for i=1:nb_sizes
    print("progress: ",i/nb_sizes*100,"%\n")
    for j=1:nb_window
        for k=1:nb_atoms
            if (j-1)*time_window < times[i] && times[i] < j*time_window && sizes[i] == k
                hist2d[j,k] += 1
            end
        end
    end
end

# Normalization
for j=1:nb_window
    count=0
    for k=1:nb_atoms
        count += hist2d[j,k]
    end
    hist2d[j,:] /= count
end

# Creating hsitory
hist_file=open(string("/home/moogmt/hist2d_size_time_",volume,"_",temperature,"_",cut_off,"_tw-",time_window,".dat"),"w")
for j=0:nb_window-1
    for k=1:nb_atoms
        write(hist_file,string(j*time_window+time_window/2," ",k," ",hist2d[j+1,k],"\n"))
    end
    write(hist_file,"\n")
end
close(hist_file)
