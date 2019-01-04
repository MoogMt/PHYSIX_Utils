GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"cubefile.jl"))
include(string(GPfolder,"clustering.jl"))

include(string(CO2folder,"markovCO2.jl"))

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

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

# distance_data=[]
# elf_data=[]
# density_data=[]
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
			if cell_mod.distance(atoms,cell,carbon1,carbon2) < 1.8
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


# function makeHistogram2D( data_x::Vector{T1}, data_y::Vector{T2}, nb_box::Vector{T3} ) where { T1 <: Real, T2<:Real, T3 <: Int }
# 	size_data
# end

# nb_box=2001
# delta_denself=1/nb_box
# min_dist=100
# max_dist=0
# for i=1:size(distance_data)[1]
#     if distance_data[i] < min_dist
#         global min_dist=distance_data[i]
#     end
#     if distance_data[i] > max_dist
#         global max_dist=distance_data[i]
#     end
# end
# delta_dist=(max_dist-min_dist)/nb_box
#
# hist2D=zeros(nb_box,nb_box)
# count_v=0
# for i=1:size(distance_data)[1]
#     for k=1:nb_box
#         if distance_data[i] > min_dist+(k-1)*delta_dist &&  distance_data[i] < min_dist+k*delta_dist
#             for l=1:nb_box
#                 if elf_data[i] > (l-1)*delta_denself && elf_data[i] < l*delta_denself
#                     hist2D[k,l] += 1
#                     break
#                 end
#             end
#             break
#         end
#     end
#     global count_v += 1
# end
#
# hist2D /= count_v
#
# file_out=open(string(folder_base,"All-elf_hist-",nb_steps,".dat"),"w")
# for i=1:nb_box
# 	for j=1:nb_box
# 		write(file_out,string(min_dist+i*delta_dist," ",j*delta_denself," ",hist2D[i,j],"\n"))
# 	end
# 	write(file_out,string("\n"))
# end
# close(file_out)
#
# hist2D=zeros(nb_box,nb_box)
# count_v=0
# for i=1:size(distance_data)[1]
#     for k=1:nb_box
#         if distance_data[i] > min_dist+(k-1)*delta_dist &&  distance_data[i] < min_dist+k*delta_dist
#             for l=1:nb_box
#                 if density_data[i] > (l-1)*delta_denself && density_data[i] < l*delta_denself
#                     hist2D[k,l] += 1
#                     break
#                 end
#             end
#             break
#         end
#     end
#     global count_v += 1
# end

# hist2D /= count_v
#
# file_out=open(string(folder_base,"All-density_hist-",nb_steps,".dat"),"w")
# for i=1:nb_box
# 	for j=1:nb_box
# 		write(file_out,string(min_dist+i*delta_dist," ",j*delta_denself," ",hist2D[i,j],"\n"))
# 	end
# 	write(file_out,string("\n"))
# end
# close(file_out)
