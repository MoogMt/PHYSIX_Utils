GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

start_=1
stride_=1
nb_steps=900

start_C=1
nbC=32
start_O=32
nbO=64

nb_atoms=nbC+nbO

nb_box=50
nb_elf=50

hist2d_out_CC=zeros(Real,nb_box,nb_elf)
hist2d_out_CC_clean=zeros(Real,nb_box,nb_elf)
hist2d_out_CC_single=zeros(Real,nb_box,nb_elf)
hist2d_out_CC_double=zeros(Real,nb_box,nb_elf)
hist2d_in_CC=zeros(Real,nb_box,nb_elf)
hist2d_out_CO=zeros(Real,nb_box,nb_elf)
hist2d_in_CO=zeros(Real,nb_box,nb_elf)

nb_points=50

cut_off_distance_clean=1.6

cut_off_distance_in=1.75
min_distance_in=1.0
delta_distance_in=(cut_off_distance_in-min_distance_in)/nb_box

cut_off_distance_in_CC=1.8
min_distance_in=1.0
delta_distance_in_CC=(cut_off_distance_in_CC-min_distance_in)/nb_box

cut_off_distance_out=2.8
min_distance_out=2.0
delta_distance_out=(cut_off_distance_out-min_distance_out)/nb_box

delta_elf=1/nb_elf

for step=start_:stride_:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,step,"_elf.cube") )
	cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
	atoms.positions=cell_mod.wrap(atoms.positions,cell)
	for atom1=start_C:start_C+nbC-1
		for atom2=start_C:start_C+nbC-1
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance < cut_off_distance_out && distance > min_distance_out  && atom1 != atom2
				distances,elfs=traceLine( atom1, atom2, nb_points, elf, atoms, cell)
				for i=1:nb_points
					ny=Int(trunc( elfs[i]/delta_elf )+1)
					hist2d_out_CC[i,ny] += 1
				end
				for atom3=start_O:start_O+nbO-1
					if cell_mod.distance(atoms.positions,cell, atom1 , atom3) > cut_off_distance_clean || cell_mod.distance(atoms.positions,cell, atom2 , atom3) > cut_off_distance_clean
						for i=1:nb_points
							ny=Int(trunc( elfs[i]/delta_elf )+1)
							hist2d_out_CC_clean[i,ny] += 1
						end
					else
						for atom4=start_O:start_O+nbO-1
							if atom3 != atom4
								if cell_mod.distance(atoms.positions,cell, atom1 , atom4) < cut_off_distance_clean && cell_mod.distance(atoms.positions,cell, atom2 , atom4) < cut_off_distance_clean
									for i=1:nb_points
										ny=Int(trunc( elfs[i]/delta_elf )+1)
										hist2d_out_CC_double[i,ny] += 1
									end
								else
									for i=1:nb_points
										ny=Int(trunc( elfs[i]/delta_elf )+1)
										hist2d_out_CC_single[i,ny] += 1
									end
								end
							end
						end
					end
				end
			elseif distance < cut_off_distance_in_CC && min_distance_in < distance && atom1 != atom2
				distances,elfs=traceLine( atom1, atom2, nb_points, elf, atoms, cell)
				for i=1:nb_points
					nx=i
					ny=Int(trunc( elfs[i]/delta_elf )+1)
					hist2d_in_CC[nx,ny] += 1
				end
			end
		end
		for atom2=start_O:start_O+nbO-1
			distance=cell_mod.distance(atoms.positions,cell, atom1 , atom2)
			if distance < cut_off_distance_out && min_distance_out < distance && atom1 != atom2
				distances,elfs=traceLine( atom1, atom2, nb_points, elf, atoms, cell)
				for i=1:nb_points
					nx=i
					ny=Int(trunc( elfs[i]/delta_elf )+1)
					hist2d_out_CO[nx,ny] += 1
				end
			elseif distance < cut_off_distance_in && min_distance_in < distance && atom1 != atom2
				distances,elfs=traceLine( atom1, atom2, nb_points, elf, atoms, cell)
				for i=1:nb_points
					nx=i
					ny=Int(trunc( elfs[i]/delta_elf )+1)
					hist2d_in_CO[nx,ny] += 1
				end
			end
		end
	end
end

for i=1:nb_points
	hist2d_out_CC[i,:] /= sum(hist2d_out_CC[i,:])
end
for i=1:nb_points
	hist2d_out_CC_clean[i,:] /= sum(hist2d_out_CC_clean[i,:])
end
for i=1:nb_points
	hist2d_out_CC_single[i,:] /= sum(hist2d_out_CC_single[i,:])
end
for i=1:nb_points
	hist2d_out_CC_double[i,:] /= sum(hist2d_out_CC_double[i,:])
end
for i=1:nb_points
	hist2d_in_CC[i,:] /= sum(hist2d_in_CC[i,:])
end
for i=1:nb_points
	hist2d_in_CO[i,:] /= sum(hist2d_in_CO[i,:])
end
for i=1:nb_points
	hist2d_out_CO[i,:] /= sum(hist2d_out_CO[i,:])
end

file_out=open(string(folder_base,"ELF_line_CC_clean.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_out_CC_clean[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CC_single.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_out_CC_single[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CC_double.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_out_CC_double[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CC_out.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_out_CC[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CC_in.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_in_CC[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CO_out.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_out_CO[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)

file_out=open(string(folder_base,"ELF_line_CO_in.dat"),"w")
for i=1:nb_elf
    for j=1:nb_points
        write(file_out,string(i/50," ",j*delta_elf," ",hist2d_in_CO[i,j],"\n"))
    end
    write(file_out,"\n")
end
close(file_out)
