GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Andrea/Glycine_TS_ELF/"

start_=1
stride_=10
nb_steps=1351

nbC=2
nbO=80
nbN=1
nbH=164


nb_atoms=nbC+nbO+nbN+nbH


file_out=open(string(folder_base,"TRAJEC_wrapped_2.xyz"),"w")
for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	if ! isfile( string(folder_base,"ELF_shoot2_",step,".cube") )
		continue
	end
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,"ELF_shoot2_",step,".cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)
	write(file_out,string(nb_atoms,"\n"))
	write(file_out,string("STEP ",step,"\n"))
	for carbon=1:nb_atoms
		if atoms.names[carbon] == "6"
			write(file_out,string("C "))
			for i=1:3
				write(file_out,string(atoms.positions[carbon,i]," "))
			end
			write(file_out,string("\n"))
		end
	end
	for oxygen=1:nb_atoms
		if atoms.names[oxygen] == "8"
			write(file_out,string("O "))
			for i=1:3
				write(file_out,string(atoms.positions[oxygen,i]," "))
			end
			write(file_out,string("\n"))
		end
	end
	for nitrogen=1:nb_atoms
		if atoms.names[nitrogen] == "7"
			write(file_out,string("N "))
			for i=1:3
				write(file_out,string(atoms.positions[nitrogen,i]," "))
			end
			write(file_out,string("\n"))
		end
	end
	for hydrogen=1:nb_atoms
		if atoms.names[hydrogen] == "1"
			write(file_out,string("H "))
			for i=1:3
				write(file_out,string(atoms.positions[hydrogen,i]," "))
			end
			write(file_out,string("\n"))
		end
	end
end
close(file_out)
