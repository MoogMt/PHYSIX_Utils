GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

start_=1
stride_=1
nb_steps=200

nb_box_distance=200
cut_off_distance=4
min_distance=0.5
delta_distance=(cut_off_distance-min_distance)/nb_box_distance

hist1d1=zeros(Real,nb_box_distance)
hist1d2=zeros(Real,nb_box_distance)
hist1d3=zeros(Real,nb_box_distance)
hist1d4=zeros(Real,nb_box_distance)

hist1dC=zeros(Real,nb_box_distance)

start_C=1
nbC=32
start_O=33
nbO=64

nb_atoms=nbC+nbO

for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)

	for atom1=start_C:start_C+nbC-1
		distances=zeros(nbO)
		for atom2=start_O:start_O+nbO-1
			distances[atom2-nbC]=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
		end
		for oxygen=1:nbO-1
			for oxygen2=oxygen+1:nbO
				if distances[oxygen] > distances[oxygen2]
					clustering.swap(distances,oxygen,oxygen2)
				end
			end
		end
		n1=Int(trunc((distances[1]-min_distance)/delta_distance)+1)
		n2=Int(trunc((distances[2]-min_distance)/delta_distance)+1)
		n3=Int(trunc((distances[3]-min_distance)/delta_distance)+1)
		n4=Int(trunc((distances[4]-min_distance)/delta_distance)+1)
		hist1d1[n1] += 1
		hist1d2[n2] += 1
		hist1d3[n3] += 1
		hist1d4[n4] += 1
	end

	for atom1=start_C:start_C+nbC-1
		distances=zeros(nbC)
		for atom2=start_C:start_C+nbC-1
			distances[atom2]=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
		end
		for carbon=1:nbC-1
			for carbon2=carbon+1:nbC
				if distances[carbon] > distances[carbon2]
					clustering.swap(distances,carbon,carbon2)
				end
			end
		end
		n1=Int(trunc((distances[2]-min_distance)/delta_distance)+1)
		hist1dC[n1] += 1
	end

end

hist1d1 /= sum(hist1d1)
hist1d2 /= sum(hist1d2)
hist1d3 /= sum(hist1d3)
hist1d4 /= sum(hist1d4)
hist1dC /= sum(hist1dC)

file_out=open(string(folder_base,"histo_base.dat"),"w")
for i=1:nb_box_distance
	write(file_out,string(i*delta_distance+min_distance," ",hist1d1[i]," ",hist1d2[i]," ",hist1d3[i]," ",hist1d4[i]," ",hist1dC[i],"\n"))
end
close(file_out)
