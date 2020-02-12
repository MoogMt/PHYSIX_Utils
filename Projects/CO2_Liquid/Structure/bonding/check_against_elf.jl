GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")
include(string(CO2folder,"markovCO2.jl"))

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1 # 0.1% cut-off for a state to be considered statististically viable

min_lag=1    # min tau
max_lag=5001 # max tau
d_lag=5      # delta tau
unit=0.005   # units of the simulation

V=8.82
nb_steps=300

n_dim=4
nb_carbon=Int(nbC*nb_steps)
distance_configurations=zeros( nb_carbon ,n_dim )
elf_configurations=zeros( nb_carbon, n_dim )
count_v=1
for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	density = cube_mod.readCube( string(folder_base,step,"_density.cube") )[3]
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms=cell_mod.wrap(atoms,cell)
	for i=1:nbC
		distances_sort = zeros( Real, nbO )
		elf_sort       = zeros( Real, nbO )
		count2=1
		for j=1:nbO
			distances_sort[count2]=cell_mod.distance(atoms,cell,i,nbC+j)
			elf_sort[count2] = cube_mod.dataInTheMiddleWME( atoms, cell , i, nbC+j, elf )
			count2+=1
		end
		for j=1:nbC-1
			for k=j+1:nbO
				if distances_sort[j] > distances_sort[k]
					stock=distances_sort[j]
					distances_sort[j]=distances_sort[k]
					distances_sort[k]=stock
					stock=elf_sort[j]
					elf_sort[j]=elf_sort[k]
					elf_sort[k]=stock
				end
			end
		end
		distance_configurations[count_v,:]=distances_sort[1:n_dim]
		elf_configurations[count_v,:]=elf_sort[1:n_dim]
		global count_v += 1
	end
end

cl_dist, icl_dist = clustering.densityPeakClusteringTrain( distance_configurations , 0.005)

size_origin=Int(nbC*nb_steps*n_dim)

max_step_train=5000
size_train=size_origin
if size_origin < max_step_train
	global size_train = max_step_train
end

elf_distance_conf=zeros(Real, size_train,2)
count_v=0
for carbon=1:n_cartbon
	for dim=1:n_dim
		elf_distance_conf[count_v,1] = distance_configuration[carbon,dim]
		elf_distance_conf[count_v,2] = elf_configuration[carbon,dim]
		global count_v += 1
	end
end
elf_configuration=[]

cl_dist, icl_dist = clustering.densityPeakClusteringTrain( elf_distance_conf , 0.005)

for i=1:size(icl_dist)[1]
	file_cluster=open(string(folder_base,"distance-cluster",i,".dat"),"w")
	for j=1:size(cl_dist)[1]
		if cl_dist[j] == i
			for k=1:n_dim
				write(file_cluster,string(distance_configurations[j,k]," "))
			end
			coord_elf=0
			for k=1:n_dim
				if elf_configurations[j,k] > 0.7 && distance_configurations[j,k] < 1.8
					coord_elf += 1
				end
			end
			write(file_cluster,string(coord_elf,"\n"))
		end
	end
	close(file_cluster)
end

cl_elf, icl_elf = clustering.densityPeakClusteringTrain( elf_configurations , 0.005)

for coord=1:n_dim
	file_elf=open(string(folder_base,"elf_vs_distance-",coord,".dat"),"w")
	for i=1:size(elf_configurations)[1]
		write(file_elf,string(i," ",distance_configurations[i,coord]," ",elf_configurations[i,coord],"\n"))
	end
	close(file_elf)
end

file_test=open(string(folder_base,"test.dat"),"w")
for i=1:size(cl)[1]
	write(file_test,string(i," ",cl_dist[i]," ",cl_elf[i],"\n"))
end
close(file_test)

for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms=cell_mod.wrap(atoms,cell)
	for carbon=1:nbC
		for oxygen=1:nbO
			if cell_mod.distance(atoms,cell,carbon1,carbon2) < 2.0
				print(step," ",carbon1," ",carbon2,"\n")
			end
		end
	end
end

nb_points=50
nb_steps=1000
elf_bond_store_double=zeros(nb_points,0)
elf_bond_store_single=zeros(nb_points,0)
elf_bond_store_interm=zeros(nb_points,0)
elf_bond_store_close=zeros(nb_points,0)
for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	for carbon=1:nbC
		for oxygen=1:nbO
			if cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 1.4
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_double = hcat( elf_bond_store_double, elfs)
			elseif cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 1.7
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_single = hcat( elf_bond_store_single, elfs)
			elseif cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 1.8
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_interm = hcat( elf_bond_store_interm, elfs)
			elseif cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 2.0
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_close = hcat( elf_bond_store_close, elfs )
			end
		end
	end
end

file_elf=open(string(folder_base,"elf_CO_double_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_double)[2]
		write(file_elf,string(elf_bond_store_double[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base,"elf_CO_single_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_single)[2]
		write(file_elf,string(elf_bond_store_single[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base,"elf_CO_interm_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_interm)[2]
		write(file_elf,string(elf_bond_store_interm[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base,"elf_CO_close_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_close)[2]
		write(file_elf,string(elf_bond_store_close[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)

nb_points=50
nb_steps=1000
elf_bond_store_bond=zeros(nb_points,0)
elf_bond_store_close1=zeros(nb_points,0)
elf_bond_store_close2=zeros(nb_points,0)
for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	for carbon1=1:nbC-1
		for carbon2=carbon1+1:nbC
			if cell_mod.distance(atoms,cell,carbon1,carbon2) < 1.75
				elfs = cube_mod.traceLine( carbon1, carbon2, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_bond = hcat( elf_bond_store_bond, elfs)
			elseif cell_mod.distance(atoms,cell,carbon1,carbon2) < 2.0
				elfs = cube_mod.traceLine( carbon1, carbon2, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_close1 = hcat( elf_bond_store_close1, elfs)
			elseif cell_mod.distance(atoms,cell,carbon1,carbon2) < 2.2
				elfs = cube_mod.traceLine( carbon1, carbon2, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_close2 = hcat( elf_bond_store_close2, elfs )
			end
		end
	end
end

file_elf=open(string(folder_base,"elf_CC_bond_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_bond)[2]
		write(file_elf,string(elf_bond_store_bond[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base,"elf_CC_close1_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_close1)[2]
		write(file_elf,string(elf_bond_store_close1[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base,"elf_CC_close2_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_close2)[2]
		write(file_elf,string(elf_bond_store_close2[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)

nb_box_line=100
delta=1/nb_box_line
hist2D=zeros(nb_points,nb_box_line)
for point=1:nb_points
	for i=1:size(elf_bond_store_bond)[2]
		for box=1:nb_box_line
			if elf_bond_store_bond[point,i] > (box-1)*delta  && elf_bond_store_bond[point,i] < box*delta
				hist2D[point,box] += 1
				break
			end
		end
	end
	hist2D[point,:]/=sum(hist2D[point,:])
end
file_hist=open(string(folder_base,"hist_CC_bond.dat"),"w")
for i=1:nb_points
	for j=1:nb_box_line
		write(file_hist,string(i/nb_points," ",j*delta," ",hist2D[i,j],"\n"))
	end
	write(file_hist,string("\n"))
end
close(file_hist)
hist2D=zeros(nb_points,nb_box_line)
for point=1:nb_points
	for i=1:size(elf_bond_store_close1)[2]
		for box=1:nb_box_line
			if elf_bond_store_close1[point,i] > (box-1)*delta  && elf_bond_store_close1[point,i] < box*delta
				hist2D[point,box] += 1
				break
			end
		end
	end
	hist2D[point,:]/=sum(hist2D[point,:])
end
file_hist=open(string(folder_base,"hist_CC_close1.dat"),"w")
for i=1:nb_points
	for j=1:nb_box_line
		write(file_hist,string(i/nb_points," ",j*delta," ",hist2D[i,j],"\n"))
	end
	write(file_hist,string("\n"))
end
close(file_hist)
hist2D=zeros(nb_points,nb_box_line)
for point=1:nb_points
	for i=1:size(elf_bond_store_close2)[2]
		for box=1:nb_box_line
			if elf_bond_store_close2[point,i] > (box-1)*delta  && elf_bond_store_close2[point,i] < box*delta
				hist2D[point,box] += 1
				break
			end
		end
	end
	hist2D[point,:]/=sum(hist2D[point,:])
end
file_hist=open(string(folder_base,"hist_CC_close2.dat"),"w")
for i=1:nb_points
	for j=1:nb_box_line
		write(file_hist,string(i/nb_points," ",j*delta," ",hist2D[i,j],"\n"))
	end
	write(file_hist,string("\n"))
end
close(file_hist)


folder_base2="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/9.8/Trajectory2/"
nb_points=50
max_step=150
elf_bond_store_bond=zeros(nb_points,0)
elf_bond_store_close=zeros(nb_points,0)
for step=1:max_step
	print("Progress: ",step/max_step*100,"%\n")
	if ! isfile(string(folder_base2,step,"_elf.cube"))
		continue
	end
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base2,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	for carbon=1:nbC
		for oxygen=1:nbO
			if cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 1.75
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_bond = hcat( elf_bond_store_bond, elfs)
			elseif cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 2.0
				elfs = cube_mod.traceLine( carbon, nbC+oxygen, nb_points , elf    , atoms , cell )[2]
				global elf_bond_store_close = hcat( elf_bond_store_close, elfs )
			end
		end
	end
end
file_elf=open(string(folder_base2,"elf_CO_bond_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_bond)[2]
		write(file_elf,string(elf_bond_store_bond[point,bond]," "))
	end
	write(file_elf,string("\n"))
end
close(file_elf)
file_elf=open(string(folder_base2,"elf_CO_close_line.dat"),"w")
for point=1:nb_points
	write(file_elf,string( (point-1)/nb_points," "))
	for bond=1:size(elf_bond_store_close)[2]
		write(file_elf,string(elf_bond_store_close[point,bond]," "))
	end
	write(file_elf,string("\n"))
end

nb_box_line=100
delta=1/nb_box_line
hist2D=zeros(nb_points,nb_box_line)
for point=1:nb_points
	for i=1:size(elf_bond_store_close1)[2]
		for box=1:nb_box_line
			if elf_bond_store_bond[point,i] > (box-1)*delta  && elf_bond_store_bond[point,i] < box*delta
				hist2D[point,box] += 1
				break
			end
		end
	end
	hist2D[point,:]/=sum(hist2D[point,:])
end
file_hist=open(string(folder_base2,"hist_CO_bond.dat"),"w")
for i=1:nb_points
	for j=1:nb_box_line
		write(file_hist,string(i/nb_points," ",j*delta," ",hist2D[i,j],"\n"))
	end
	write(file_hist,string("\n"))
end
close(file_hist)
hist2D=zeros(nb_points,nb_box_line)
for point=1:nb_points
	for i=1:size(elf_bond_store_close)[2]
		for box=1:nb_box_line
			if elf_bond_store_close[point,i] > (box-1)*delta  && elf_bond_store_close[point,i] < box*delta
				hist2D[point,box] += 1
				break
			end
		end
	end
	hist2D[point,:]/=sum(hist2D[point,:])
end
file_hist=open(string(folder_base2,"hist_CO_close.dat"),"w")
for i=1:nb_points
	for j=1:nb_box_line
		write(file_hist,string(i/nb_points," ",j*delta," ",hist2D[i,j],"\n"))
	end
	write(file_hist,string("\n"))
end
close(file_hist)


elfs_data=zeros(Real,max_step,2)
file_out=open(string(folder_base2,"elf_vs_distance.dat"),"w")
count_v=1
for step=1:max_step
	print("Progress: ",step/max_step*100,"%\n")
	if ! isfile(string(folder_base2,step,"_elf.cube"))
		continue
	end
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base2,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	for carbon=1:nbC
		for oxygen=1:nbO
			if cell_mod.distance(atoms,cell,carbon,nbC+oxygen) < 2.5
				elfs = cube_mod.dataInTheMiddleWME( atoms, cell , carbon , nbC+oxygen, elf )
				write(file_out,string( cell_mod.distance(atoms,cell,carbon,nbC+oxygen), " ", elfs, "\n" ) )
				if  count_v <= max_step_train
					elfs_data[count_v,1] = elfs
					elfs_data[count_v,2] = cell_mod.distance(atoms,cell,carbon,nbC+oxygen)
					global count_v += 1
				end
			end
		end
	end
end
close(file_out)

max_dist = 0
min_dist = 10000
for i=1:max_step
	if max_dist < elfs_data[i,2]
		global max_dist = elfs_data[i,2]
	end
	if min_dist > elfs_data[i,2]
		global min_dist = elfs_data[i,2]
	end
end

nb_box_dist=100
delta_distance=(max_dist-min_dist)/nb_box_dist

nb_box_elf=100
delta_elf=1/nb_box_elf

hist2D=zeros(nb_box_dist,nb_box_elf)
for k=1:size(elfs_data)[1]
	for i=1:nb_box_dist
		if elfs_data[k,2] > min_dist+(i-1)*delta_distance && elfs_data[k,2] < min_dist+i*delta_distance
			for j=1:nb_box_elf
				if elfs_data[k,1] > (j-1)*delta_elf && elfs_data[k,1] < j*delta_elf
					hist2D[i,j] += 1
				end
			end
		end
	end
end
hist2D /= sum(hist2D)

file_out=open(string(folder_base2,"hist_elf_vs_dist.dat"),"w")
for i=1:nb_box_dist
	for j=1:nb_box_elf
		write(file_out,string(i*delta_distance," ",j*delta_elf," ",hist2D[i,j],"\n"))
	end
	write(file_out,"\n")
end
close(file_out)


# Clustering
max_step_train=4000
elfs_data=elfs_data[1:max_step_train]
cl, icl = clustering.densityPeakClusteringTrain( elfs_data , 0.01)
for i=1:size(icl)[1]
	file_out=open(string(folder_base2,"elf_clusterCO-",i,".dat"),"w")
	for j=1:size(cl)[1]
		if cl[j] == i
			write(file_out,string(elfs_data[j,2]," ",elfs_data[j,1],"\n"))
		end
	end
	close(file_out)
end
