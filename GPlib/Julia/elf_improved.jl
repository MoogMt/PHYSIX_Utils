include("atoms.jl")
include("cell.jl")
include("cubefile.jl")


unit=0.005
stride=5


for step=1:10:1001

atoms, cell1, ELF1 = cube_mod.readCube(string("/media/moogmt/Stock/cube_TS14_gly/ELF_shoot1_",step,".cube"))

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

atom1s=[]
atom2s=[]
file_in=open("/home/moogmt/log_weird.dat")
lines=readlines(file_in);
close(file_in)
for i=1:size(lines)[1]
    line=split(lines[i])
    step=parse(Int64,line[1])
    if step == 1
        push!(atom1s,parse(Int64,line[2]))
        push!(atom2s,parse(Int64,line[3]))
    end
end

n_dots=10

# Aligning ELF
for atom=1:size(atoms.names)[1]
	for i=1:3
		atoms.positions[atom,i]=cell_mod.wrap(atoms.positions[atom,i],params[i])
	end
	atoms.positions[atom,:]-=ELF1.origin
end


file=open(string("/home/moogmt/test.xyz"),"w")
file_out=open(string("/home/moogmt/test_dist.dat"),"w")
write(file,string(2*size(atom1s)[1],"\nCHECK\n"))
for w=1:size(atom1s)[1]
atom1=atom1s[w]
atom2=atom2s[w]

distanceatm=cell_mod.distance( atoms, cell2, atom1, atom2 )
# Ajusting atom2
for j=1:3
    di = atoms.positions[atom1,j]-atoms.positions[atom2,j]
    if di > cell2.length[j]*0.5
	atoms.positions[atom2,j] += cell2.length[j]
    end
    if di < -cell2.length[j]*0.5
	atoms.positions[atom2,j] -= cell2.length[j]
    end
end
temp1=atoms.positions[atom1,:]
temp2=atoms.positions[atom1,:]
for i=1:3
	temp1[i]=cell_mod.wrap(atoms.positions[atom1,i],params[i])
	temp2[i]=cell_mod.wrap(atoms.positions[atom2,i],params[i])
end
center=(atoms.positions[atom1,:]+atoms.positions[atom2,:])/2.
# Wrapping PBC
for i=1:3
	center[i]=cell_mod.wrap(center[i], params[i])
end

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

distance1=0
for i=1:3
    distance1+=cell_mod.dist1D( index[i]*cell2.length[i]/ELF1.nb_vox[i], center[i], cell2.length[i] )^2
end
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

for i=1:3
	index[i] += 1
end

distance3=0
for i=1:3
	distance3 += cell_mod.dist1D(temp1[i],index[i]*params[i]/ELF1.nb_vox[i],params[i])^2
end
distance3 = sqrt(distance3)
if distance3 > distanceatm
	print(atom1," ",atom2,"\n")
end


write(file_out,string(distance3," ",ELF1.matrix[index[1],index[2],index[3]],"\n"))

end
end

close(file)
close(file_out)
