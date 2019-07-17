GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

start_=1
stride_=1
nb_steps=900

nb_box_distance=100
cut_off_distance=6
min_distance=0.5
delta_distance=(cut_off_distance-min_distance)/nb_box_distance

nb_box_elf=50
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
		for atom2=start_C:start_C+nbC
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

file_out=open(string(folder_base,"histCC.dat"),"w")
for i=1:nb_box_distance
	for j=1:nb_box_elf
		write(file_out,string(i*delta_distance+min_distance*0.529177," ",j*delta_elf," ",hist2d[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)

min_distance_test=2 #Bohr
cut_off_distance_test=6
nb_steps=100

nb_box_distance=100
cut_off_distance=6
delta_distance=(cut_off_distance_test-min_distance_test)/nb_box_distance

nb_box_elf=50
delta_elf=1/nb_box_elf

hist1d=zeros(Int,nb_box_distance,nb_box_elf)

for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_C:start_C+nbC
		for atom2=start_C:start_C+nbC
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance > min_distance_test && cut_off_distance_test > distance && atom1 != atom2
				elf_value = cube_mod.dataInTheMiddleWME( atoms, cell , atom1, atom2, elf )
				if elf_value > 0.5
					dist=8.82/0.529177
					for oxygen=start_O:start_O+nbO
						dist_local = cell_mod.distance()
						if dist > dist_local
							dist =  dist_local
						end
					end
					hist1d[ Int(trunc()) ] += 1
				end
			end
		end
	end

end


file_out=open(string(folder_base,"histCC_confirmation.dat"),"w")
for i=1:nb_box_distance
	for j=1:nb_box_elf
		write(file_out,string((i*delta_distance+min_distance_test+min_distance_test)*0.529177," ",j*delta_elf," ",hist2d[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)
