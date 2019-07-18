GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

start_=1
stride_=1
nb_steps=950

nb_box_distance=100
cut_off_distance=2.5
min_distance=0.7
delta_distance=(cut_off_distance-min_distance)/nb_box_distance

nb_box_elf=100
delta_elf=1/nb_box_elf

hist2d=zeros(Int,nb_box_distance,nb_box_elf)

start_C=1
nbC=32
start_O=32
nbO=64

nb_atoms=nbC+nbO

for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_C:start_C+nbC
		for atom2=start_O:start_O+nbO
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance < cut_off_distance && min_distance < distance && atom1 != atom2
				elf_value = cube_mod.dataInTheMiddleWME( atoms, cell , atom1, atom2, elf )
				nx=Int(trunc((distance-min_distance)/delta_distance))+1
				ny=Int(trunc(elf_value/delta_elf))+1
				hist2d[nx,ny] += 1
			end
		end
	end

end

file_out=open(string(folder_base,"histCO.dat"),"w")
for i=1:nb_box_distance
	for j=1:nb_box_elf
		write(file_out,string(i*delta_distance+min_distance," ",j*delta_elf," ",hist2d[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)

#---------------------------------------------------------------------------

start_=1
stride_=1
nb_steps=100

neighbor=4


start_C=1
nbC=32
start_O=33
nbO=64


file_out=open(string(folder_base,"1-4_ELF.dat"),"w")
file_out2=open(string(folder_base,"1-4_distance.dat"),"w")
for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_C:start_C+nbC-1
		index_=clustering.simpleSequence(nbO).+nbC
		distances=zeros(nbO)
		for atom2=start_O:start_O+nbO-1
			distances[atom2-nbC]=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
		end
		for i=1:nbO
			for j=i+1:nbO
				if distances[i] > distances[j]
					clustering.swap(distances,i,j)
					clustering.swap(index_,i,j)
				end
			end
		end
		write(file_out,string(step," "))
		write(file_out2,string(step," "))
		for i=1:neighbor
			write(file_out,string(cube_mod.dataInTheMiddleWME( atoms, cell , atom1, index_[i], elf )," "))
			write(file_out2,string(distances[i]," "))
		end
		write(file_out,string("\n"))
		write(file_out2,string("\n"))
	end
end
close(file_out)
close(file_out2)

#---------------------------------------------------------------------------

start_=1
stride_=1
nb_steps=100

nb_box_distance=50
cut_off_distance=2.5 #bohr?
min_distance=0 #Bohr?
cut_histo=8
delta_distance=(cut_histo)/nb_box_distance

start_C=1
nbC=32
start_O=32
nbO=64

hist2d_1=zeros(Int,nb_box_distance,nb_box_distance)
hist2d_2=zeros(Int,nb_box_distance,nb_box_distance)

for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_C:start_C+nbC
		for atom2=start_C:start_C+nbC
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance > min_distance && cut_off_distance > distance && atom1 != atom2
				elf_value = cube_mod.dataInTheMiddleWME( atoms, cell , atom1, atom2, elf )
				if elf_value > 0.5
					dist=8.82
					target=0
					for oxygen=start_O:start_O+nbO
						middle=(atoms.positions[atom1,:].+atoms.positions[atom2,:]).*0.5
						dist_local = cell_mod.distance(atoms.positions[oxygen,:],middle,cell.length[:])
						if dist > dist_local
							dist =  dist_local
							target=oxygen
						end
					end
					nx=Int(trunc(cell_mod.distance(atoms.positions,cell,atom1,target)/delta_distance))+1
					ny=Int(trunc(cell_mod.distance(atoms.positions,cell,atom2,target)/delta_distance))+1
					hist2d_1[ nx, ny ] += 1
				else
					dist=8.82
					target=0
					for oxygen=start_O:start_O+nbO
						middle=(atoms.positions[atom1,:].+atoms.positions[atom2,:]).*0.5
						dist_local = cell_mod.distance(atoms.positions[oxygen,:],middle,cell.length[:])
						if dist > dist_local
							dist =  dist_local
							target=oxygen
						end
					end
					nx=Int(trunc(cell_mod.distance(atoms.positions,cell,atom1,target)/delta_distance))+1
					ny=Int(trunc(cell_mod.distance(atoms.positions,cell,atom2,target)/delta_distance))+1
					hist2d_2[ nx, ny ] += 1
				end
			end
		end
	end

end


file_out_1=open(string(folder_base,"histCC_DistanceO_1.dat"),"w")
file_out_2=open(string(folder_base,"histCC_DistanceO_2.dat"),"w")
for i=1:nb_box_distance
	for j=1:nb_box_distance
		write(file_out_1,string(i*delta_distance," ",j*delta_distance," ",hist2d_1[i,j],"\n"))
		write(file_out_2,string(i*delta_distance," ",j*delta_distance," ",hist2d_2[i,j],"\n"))
	end
	write(file_out_1,string("\n"))
	write(file_out_2,string("\n"))
end
close(file_out_1)
close(file_out_2)
