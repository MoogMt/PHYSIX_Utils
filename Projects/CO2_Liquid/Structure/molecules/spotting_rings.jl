GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz


V=9.5
T=3000
folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
#folder_base="/home/moogmt/CO2/CO2_AIMD/8.82/3000K/"
file=string(folder_base,"TRAJEC_wrapped.xyz")

traj = filexyz.readFastFile( file )
cell=cell_mod.Cell_param(V,V,V)

cut_off=1.75
nbC=32
nbO=64

nb_atoms=nbC+nbO
nb_steps=size(traj)[1]

matrix_full=zeros(nb_atoms,nb_atoms)
used=zeros(nb_atoms)

file_out=open(string(folder_base,"/Data/log_rings.dat"),"w")
for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	matrix_full=zeros(nb_atoms,nb_atoms)
	# Computing matrix
	for atom1=1:nb_atoms-1
		for atom2=1:nb_atoms
			if cell_mod.distance(traj[step],cell,atom1,atom2) < cut_off
				matrix_full[atom1,atom2] = 1
				matrix_full[atom2,atom1] = 1
			end
		end
	end
	used=zeros(nb_atoms)
	# Testing rings
	for carbon=1:nbC
		if used[carbon] == 1
			continue
		end
		for carbon2=carbon+1:nbC
			if used[carbon2] == 1
				continue
			end
			for oxygen=1:nbO
				if used[oxygen] == 1
					continue
				end
				if matrix_full[carbon,nbC+oxygen] == 1 && matrix_full[carbon2,nbC+oxygen] == 1
					for carbon3=carbon2+1:nbC
						if used[carbon3] == 1
							continue
						end
						for oxygen2=1:nbO
							if used[oxygen2] == 1
								continue
							end
							if matrix_full[carbon,nbC+oxygen2] == 1 && matrix_full[carbon3,nbC+oxygen2] == 1
								for oxygen3=1:nbO
									if used[oxygen3] == 1
										continue
									end
									if matrix_full[carbon2,oxygen3] == 1 && matrix_full[carbon3,oxygen3] == 1
										Base.write(file_out,string(step," ",carbon," ",carbon2," ",carbon3," ",oxygen," ",oxygen2," ",oxygen3,"\n"))
										used[carbon]=1
										used[carbon2]=1
										used[carbon3]=1
										used[oxygen]=1
										used[oxygen2]=1
										used[oxygen3]=1
									end
								end
							end
						end
					end
				end
			end
		end
	end
end
close(file_out)
