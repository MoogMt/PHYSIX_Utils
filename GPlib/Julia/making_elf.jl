include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

atom1=28
atom2=82

file=open(string("/home/moogmt/check2.dat"),"w")

unit=0.005
stride=5

current_volume=8.82
cell2=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nelf=50
ndist=50
min_elf=0
max_elf=1
delf=(max_elf-min_elf)/nelf
min_dist=1
max_dist=3
ddist=(max_dist-min_dist)/ndist
hist2d=zeros(50,50)

count2=0
for i=0:49
    print(string(i,"\n"))
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/",i,"_structure/")
    atoms1, cell1, ELF1 = cube_mod.readCube(string(folder,"ELF.cube"))
    for atom1=1:32
    for atom2=33:96
    distanceatm=cell_mod.distance(atoms1,cell2,atom1,atom2)
    if distanceatm < 3
    temp2=atoms1.positions[atom2,:]
            # PBC for orthorombic structure
    for j=1:3
        di = atoms1.positions[atom1,j]-atoms1.positions[atom2,j]
        if di > cell2.length[j]*0.5
            atoms1.positions[atom2,j] += cell2.length[j]
        end
        if di < -cell2.length[j]*0.5
            atoms1.positions[atom2,j] -= cell2.length[j]
        end
    end
    #pos1=cube_mod.getClosest(( atoms1.positions[atom1,:]+atoms1.positions[atom2,:])/2 ,ELF1)
    # Wrapping for PBC
    position_origin=(atoms1.positions[atom1,:]+atoms1.positions[atom2,:]-ELF1.origin)/2
    for j=1:3
        position_origin[j] = cell_mod.wrap( position_origin[j] , cell2.length[j] )
    end
    # original index
    floats=[0,0,0]
    for j=1:3
        check=position_origin[j]*ELF1.nb_vox[j]/cell2.length[j]
        floats[j]=trunc(check)
        if check - floats[j] > 0.5
            floats[j] += 1
        end
        if floats[j] == 0
            floats[j]=1
        end
    end
    ddmax=0
    for j=1:3
        ddmax += ELF1.vox_vec[j,j]
    end
    elf_value=ELF1.matrix[ floats[1], floats[2], floats[3]]
    count=1
    ddmax = (ddmax/3)^2
    for p=-1:2:1
        for j=-1:2:1
            for k=-1:2:1
                new_floats=[0,0,0]
                new_floats[1]=floats[1]+p
                new_floats[2]=floats[2]+j
                new_floats[3]=floats[3]+k
                for l=1:3
                    if new_floats[l] <= 0
                        new_floats[l] = ELF1.nb_vox[l] - new_floats[l]
                    end
                    if new_floats[l] > ELF1.nb_vox[l]
                        new_floats[l] = new_floats[l]- ELF1.nb_vox[l]
                    end
                end
                distance2=0
                for l=1:3
                    distance2 += (new_floats[l]*cell2.length[l]/ELF1.nb_vox[l]+ELF1.origin[l]- position_origin[l])^2
                end
                #if distance2 < ddmax
                    elf_value += ELF1.matrix[ new_floats[1], new_floats[2], new_floats[3] ]
                    count +=1
                #end
            end
        end
    end

    elf_value = elf_value/count
    if elf_value < max_elf && distanceatm < max_dist
    for j=1:ndist
        for k=1:nelf
            if distanceatm > (j-1)*ddist + min_dist && distanceatm < j*ddist + min_dist
                if elf_value > (k-1)*delf + min_elf && elf_value < k*delf + min_elf
                    hist2d[j,k] += 1
                    count2 += 1
                end
            end
        end
    end
    end
    write(file,string( cell_mod.distance(atoms1,cell2,atom1,atom2), " ", elf_value ," ",i*unit*stride," ",i,"\n") )
    end
    end
    end
end

close(file)

filehist=open("/home/moogmt/hist_elfdist2.dat","w")
for i=1:nelf
    for j=1:ndist
        write(filehist,string( i*ddist+min_dist, " ", j*delf+min_elf," ",hist2d[i,j]/count2,"\n"))
    end
    write(filehist,"\n")
end
close(filehist)
