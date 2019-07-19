GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz


V=8.82
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"
folder_base="/home/moogmt/CO2/CO2_AIMD/8.82/3000K/"
file=string(folder_base,"TRAJEC_wrapped.xyz")

traj = filexyz.readFastFile( file )
cell=cell_mod.Cell_param(V,V,V)

nb_box_distance=100

cut_off_distance_C=3.2
min_distance_C=1.0
delta_distance_C=(cut_off_distance_C-min_distance_C)/nb_box_distance

cut_off_distance_O=2.8
min_distance_O=1.0
delta_distance_O=(cut_off_distance_O-min_distance_O)/nb_box_distance

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
nb_steps=size(traj)[1]

for step=1:nb_steps
	print("Progress: ",step/nb_steps*100,"%\n")
	for atom1=start_C:start_C+nbC-1
		distances=zeros(nbO)
		for atom2=start_O:start_O+nbO-1
			distances[atom2-nbC]=cell_mod.distance(traj[step].positions,cell, atom1 , atom2)
		end
		for oxygen=1:nbO-1
			for oxygen2=oxygen+1:nbO
				if distances[oxygen] > distances[oxygen2]
					clustering.swap(distances,oxygen,oxygen2)
				end
			end
		end
		n1=Int(trunc((distances[1]-min_distance_O)/delta_distance_O)+1)
		n2=Int(trunc((distances[2]-min_distance_O)/delta_distance_O)+1)
		n3=Int(trunc((distances[3]-min_distance_O)/delta_distance_O)+1)
		n4=Int(trunc((distances[4]-min_distance_O)/delta_distance_O)+1)
		if n1 > 0 && n1 < nb_box_distance
			hist1d1[n1] += 1
		end
		if n2 > 0 && n2 < nb_box_distance
			hist1d2[n2] += 1
		end
		if n3 > 0 && n3 < nb_box_distance
			hist1d3[n3] += 1
		end
		if n4 > 0 && n4 < nb_box_distance
			hist1d4[n4] += 1
		end
	end
	for atom1=start_C:start_C+nbC-1
		distances=zeros(nbC)
		for atom2=start_C:start_C+nbC-1
			distances[atom2]=cell_mod.distance(traj[step].positions,cell, atom1 , atom2)
		end
		for carbon=1:nbC-1
			for carbon2=carbon+1:nbC
				if distances[carbon] > distances[carbon2]
					clustering.swap(distances,carbon,carbon2)
				end
			end
		end
		n1=Int(trunc((distances[2]-min_distance_C)/delta_distance_C)+1)
		if n1 > 0 && n1 < nb_box_distance
			hist1dC[n1] += 1
		end
	end
end

hist1d1 /= sum(hist1d1)
hist1d2 /= sum(hist1d2)
hist1d3 /= sum(hist1d3)
hist1d4 /= sum(hist1d4)
hist1dC /= sum(hist1dC)

file_out=open(string(folder_base,"histo_base.dat"),"w")
for i=1:nb_box_distance
	Base.write(file_out,string(i*delta_distance+min_distance," ",hist1d1[i]," ",hist1d2[i]," ",hist1d3[i]," ",hist1d4[i]," ",hist1dC[i],"\n"))
end
close(file_out)
