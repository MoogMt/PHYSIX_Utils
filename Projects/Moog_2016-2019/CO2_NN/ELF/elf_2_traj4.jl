GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/9.0/Trajectory/"

start_=1
stride_=1
nb_steps=99

start_C=1
nbC=32
start_O=32
nbO=64

nb_atoms=nbC+nbO

file_out=open(string(folder_base,"TRAJEC_wrapped.xyz"),"w")
for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_structure/ELF.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)
	write(file_out,string(nb_atoms,"\n"))
	write(file_out,string("STEP ",step,"\n"))
	for carbon=1:nbC
		write(file_out,string("C "))
		for i=1:3
			write(file_out,string(atoms.positions[carbon,i]," "))
		end
		write(file_out,string("\n"))
	end
	for oxygen=nbC+1:nbO+nbC
		write(file_out,string("O "))
		for i=1:3
			write(file_out,string(atoms.positions[oxygen,i]," "))
		end
		write(file_out,string("\n"))
	end
end
close(file_out)
