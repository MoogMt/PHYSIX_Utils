GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"cubefile.jl"))
include(string(CO2folder,"markovCO2.jl"))

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1 # 0.1% cut-off for a state to be considered statististically viable

min_lag=1    # min tau
max_lag=5001 # max tau
d_lag=5      # delta tau
unit=0.005   # units of the simulation

V=8.82
nb_steps=1000

atom1_data=[]
atom2_data=[]
distance_data=[]
elf_data=[]
for i=1:nb_steps
    print("Progress: ",i/nb_steps*100,"%\n")
    atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,i,"_elf.cube") )
    density = cube_mod.readCube( string(folder_base,i,"_density.cube") )[3]
    cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
    atoms=cell_mod.wrap(atoms,cell)
    for i=1:nbC
        for j=nbC+1:nbC+nbO
            if cell_mod.distance(atoms,cell,i,j) < 2.5
                push!(distance_data,cell_mod.distance(atoms,cell,i,j))
                push!(elf_data,cube_mod.dataInTheMiddleWME( atoms, cell , i, j, elf ))
            end
        end
    end
end

nb_box=200
delta_denself=1/nb_box
min_dist=100
max_dist=0
for i=1:size(distance_data)[1]
    if distance_data[i] < min_dist
        global min_dist=distance_data[i]
    end
    if distance_data[i] > max_dist
        global max_dist=distance_data[i]
    end
end
delta_dist=(max_dist-min_dist)/nb_box

hist2D=zeros(nb_box,nb_box)
count_v=0
for i=1:size(distance_data)[1]
    for k=1:nb_box
        if distance_data[i] > min_dist+(k-1)*delta_dist &&  distance_data[i] < min_dist+k*delta_dist
            for l=1:nb_box
                if elf_data[i] > (l-1)*delta_denself && elf_data[i] < l*delta_denself
                    hist2D[k,l] += 1
                    break
                end
            end
            break
        end
    end
    global count_v += 1
end

hist2D /= count_v

file_out=open(string(folder_base,"elf_hist.dat"),"w")
for i=1:nb_box
	for j=1:nb_box
		write(file_out,string(min_dist+i*delta_dist," ",j*delta_denself," ",hist2D[i,j],"\n"))
	end
	write(file_out,string("\n"))
end
close(file_out)
