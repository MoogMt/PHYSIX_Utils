GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz


# Looking for C-C bonds and computing distances between C-C (exclude C-O-C triangles)

V=9.5
T=3000

cut_off=1.8
cut_off2=1.9
nbC=32
nbO=64

nb_atoms=nbC+nbO

folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder_base,"TRAJEC_wrapped.xyz")

traj = filexyz.readFastFile( file )
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]

file_out=open(string(folder_base,"Data/Trajec_CC.xyz"),"w")
file_out2=open(string(folder_base,"Data/Log_CC.xyz"),"w")
for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	count_=0
	for carbon=1:nbC
		for carbon2=carbon+1:nbC
			if cell_mod.distance(traj[step],cell,carbon,carbon2) < cut_off
				found=false
				for oxygen=1:nbO
					if cell_mod.distance(traj[step],cell,carbon,oxygen) < cut_off && cell_mod.distance(traj[step],cell,carbon2,oxygen) < cut_off
						found = true
					end
				end
				if ! false
					count_ += 1
					Base.write(file_out2,string(step," ",carbon," ",carbon2,"\n"))
				end
			end
		end
	end
	if count_ > 0
		write(file_out,string("96\n"))
		write(file_out,string("STEP\n"))
		for carbon=1:nbC
			write(file_out,string("C "))
			for k=1:3
				write(file_out,string(traj[step].positions[carbon,k]," "))
			end
			write(file_out,string("\n"))
		end
		for oxygen=nbC+1:nb_atoms
			write(file_out,string("O "))
			for k=1:3
				write(file_out,string(traj[step].positions[oxygen,k]," "))
			end
			write(file_out,string("\n"))
		end
	end
end
close(file_out)
close(file_out2)
