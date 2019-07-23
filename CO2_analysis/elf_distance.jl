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

start_C=1
nbC=32
start_O=32
nbO=64

nb_atoms=nbC+nbO

nb_box=100
nb_elf=50

hist2D_=zeros(Real,nb_box,nb_elf)

cut_off_distance=2.8
min_distance=0.8

delta_distance=(cut_off_distance-min_distance)/nb_box
delta_elf=1/nb_elf

for step=start_:stride_:nb_steps
    print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)
	for carbon=start_C:start_C+nbC-1
		for oxygen=start_O:start_O+nbO-1
			distance=cell_mod.distance( atoms.positions, cell, carbon , oxygen )
			if distance > cut_off_distance || distance < min_distance
				continue
			end
			elf_value=dataInTheMiddleWME( atoms, cell , carbon, oxygen, elf )
			nx=Int(trunc(elf_value/delta_elf)+1)
			ny=Int(trunc((distance-min_distance)/delta_distance)+1)
			hist2D_[ny,nx] += 1
		end
	end
end

sum_=sum(hist2D_[:,:])
for i=1:nb_elf
	for j=1:nb_box
		hist2D_[i,j] /= sum_
	end
end

file_out=open(string(folder_base,"elf_distance-",nb_steps,".dat"),"w")
for i=1:nb_box
	for j=1:nb_elf
		write(file_out,string(i*delta_distance+min_distance," ",j*delta_elf," ",hist2D_[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)
