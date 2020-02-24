include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

step_max=1000
step_min=1
d_step=1

distance_data=[]
elf_data=[]
used=[]

for step=step_min:d_step:step_max

#------------------------------------------------------------------------------
atoms, cell1, ELF1 = cube_mod.readCube(string("/home/moogmt/CO2/CO2_AIMD/ELF/ELF_8.82_results/",step,"_elf.cube"))
#------------------------------------------------------------------------------

#---------------
# Compute cell
#------------------------------------------------------------
params=zeros(3)
for j=1:3
    for k=1:3
        params[j] += cell1[j,k]^2
    end
    params[j] = sqrt(params[j])
end
cell2=cell_mod.Cell_param(params[1],params[2],params[3])
#------------------------------------------------------------

#----------------
# Aligning ELF
#-------------------------------------------------------------------------------
for atom=1:size(atoms.names)[1]
	for i=1:3
		atoms.positions[atom,i]=cell_mod.wrap(atoms.positions[atom,i],params[i])
	end
	atoms.positions[atom,:]-=ELF1.origin
end
#-------------------------------------------------------------------------------

file_out=open(string("/home/moogmt/test_dist.dat"),"w")

for atom1=1:32

distances_check=zeros(64)
indexes=zeros(64)

for atom3=33:96

#------------------------------------------------------------------------------
indexes[atom3-32] = atom3-32
distances_check[atom3-32] = cell_mod.distance( atoms, cell2, atom1, atom3 )
#------------------------------------------------------------------------------
# end atom2
end

for i=1:64
    for j=i+1:64
        if distances_check[j] < distances_check[i]
            stock=indexes[j]
            indexes[j] = indexes[i]
            indexes[i] = stock
        end
    end
end

atom2=Int(trunc(indexes[2]))

#-----------------
# Ajusting atom2
#------------------------------------------------------------------------------
for j=1:3
    di = atoms.positions[atom1,j]-atoms.positions[atom2,j]
    if di > cell2.length[j]*0.5
	atoms.positions[atom2,j] += cell2.length[j]
    end
    if di < -cell2.length[j]*0.5
	atoms.positions[atom2,j] -= cell2.length[j]
    end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
temp1=atoms.positions[atom1,:]
temp2=atoms.positions[atom1,:]
for i=1:3
	temp1[i]=cell_mod.wrap(atoms.positions[atom1,i],params[i])
	temp2[i]=cell_mod.wrap(atoms.positions[atom2,i],params[i])
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
center=(atoms.positions[atom1,:]+atoms.positions[atom2,:])/2.
#------------------------------------------------------------------------------
#--------------
# Wrapping PBC
#------------------------------------------------------------------------------
for i=1:3
	center[i]=cell_mod.wrap(center[i], params[i])
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
index=[0,0,0]
for i=1:3
	check=center[i]*ELF1.nb_vox[i]/params[i]
    	index[i]=trunc(check)
	if check - index[i] > 0.5
		index[i] += 1
	end
	if index[i] > ELF1.nb_vox[i]-1
		index[i]=0
	end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
distance1=0
for i=1:3
    distance1+=cell_mod.dist1D( index[i]*cell2.length[i]/ELF1.nb_vox[i], center[i], cell2.length[i] )^2
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
for l=-1:1:1
	for m=-1:1:1
		for n=-1:1:1
			new_index=[0,0,0]
			new_index[1]=index[1]+l
			new_index[2]=index[2]+m
			new_index[3]=index[3]+n
			for l=1:3
		    		if new_index[l] < 0
					new_index[l] = ELF1.nb_vox[l] - new_index[l]
		    		end
		    		if new_index[l] >= ELF1.nb_vox[l]-1
					new_index[l] = new_index[l]-ELF1.nb_vox[l]
		    		end
			end
			distance2=0
			for i=1:3
		    		distance2+=cell_mod.dist1D( new_index[i]*cell2.length[i]/ELF1.nb_vox[i], center[i], params[i] )^2
			end
			if distance2 < distance1
		    		distance1 = distance2
		    		index=new_index
			end
	    	end
	end
end
#------------------------------------------------------------------------------

#------------------
for i=1:3
	index[i] += 1
end
#------------------

#-------------
# Adding data
#---------------------------------------------------
push!( used, 0)
push!( distance_data,  cell_mod.distance( atoms, cell2, atom1, atom2 ) )
push!( elf_data, ELF1.matrix[ index[1], index[2], index[3]] )
#---------------------------------------------------

# end atom1
end

print("step:",step,"\n")
# end step
end

min_distance=0.9
max_distance=3.0
n_distance=100
d_distance=(max_distance-min_distance)/n_distance

min_elf=0
max_elf=1.0
n_elf=100
d_elf=(max_elf-min_elf)/n_elf

hist_2D=zeros(n_elf,n_distance)
count_=0
for i=1:n_elf
	print("progress: ",i/n_elf*100,"%\n")
	for j=1:n_distance
		for k=1:size(used)[1]
			if (i-1)*d_elf+min_elf < elf_data[k] && i*d_elf+min_elf > elf_data[k] && (j-1)*d_distance+min_distance < distance_data[k] && j*d_distance+min_distance > distance_data[k]
				hist_2D[i,j] += 1
				global count_ += 1
			end
		end
	end
end

hist_2D /= count_

file_hist=open("/home/moogmt/hist.dat","w")
for i=1:n_elf
    for j=1:n_distance
        write(file_hist,string(j*d_distance+min_distance," ",i*d_elf+min_elf," ",hist_2D[i,j],"\n"))
    end
    write(file_hist,"\n")
end
close(file_hist)
