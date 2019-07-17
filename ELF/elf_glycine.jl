GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod

folder_base="/media/moogmt/Stock/Andrea/Glycine_TS_ELF/"

nb_steps=1091


nb_box_distance=50
cut_off_distance=5
min_distance=0.2
delta_distance=(cut_off_distance-min_distance)/nb_box_distance

nb_box_elf=20
delta_elf=1/nb_box_elf

hist2d=zeros(Int,nb_box_distance,nb_box_elf)


start_C=1
nbC=2
start_O=3
nbO=80
start_H=5
nbH=164
start_N=84
nbN=1

nb_atoms=nbC+nbO+nbN+nbH

for step=1:10:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,"ELF_shoot1_",step,".cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_H:start_H+nbH
		for atom2=start_N:start_N+nbN
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance < cut_off_distance && min_distance < distance
				elf_value = cube_mod.dataInTheMiddleWME( atoms, cell , atom1, atom2, elf )
				nx=Int(trunc((distance-min_distance)/delta_distance))+1
				ny=Int(trunc(elf_value/delta_elf))+1
				hist2d[nx,ny] += 1
			end
		end
	end

end

file_out=open(string(folder_base,"histNH.dat"),"w")
for i=1:nb_box_distance
	for j=1:nb_box_elf
		write(file_out,string(i*delta_distance+min_distance," ",j*delta_elf," ",hist2d[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)
