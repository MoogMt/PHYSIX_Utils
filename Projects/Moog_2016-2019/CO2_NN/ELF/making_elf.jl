include("atoms.jl")
include("cell.jl")
include("cubefile.jl")


unit=0.005
stride=5

nelf=200
ndist=200
min_elf=0
max_elf=1
delf=(max_elf-min_elf)/nelf
min_dist=0.5
max_dist=3
ddist=(max_dist-min_dist)/ndist
hist2d=zeros(nelf,ndist)

file_check=open(string("/home/moogmt/error_distance.dat"),"w")
log=open(string("/home/moogmt/log_weird.dat"),"w")

count2=0
for i=1:10:1001
    print(string(i,"\n"))
    folder=string()
    atoms1, cell1, ELF1 = cube_mod.readCube(string("/media/moogmt/Stock/cube_TS14_gly/ELF_shoot1_",i,".cube"))
    #--------------------------------------
    # Compute cell parameters
    #--------------------------------------
    params=zeros(3)
    for j=1:3
        for k=1:3
            params[j] += cell1[j,k]^2
        end
        params[j] = sqrt(params[j])
    end
    cell2=cell_mod.Cell_param(params[1],params[2],params[3])
    nb_atoms=size(atoms1.names)[1]
    #--------------------------------------
    # Loop over atoms
    #-----------------------------------
    for atom1=3:84
        for atom2=85:nb_atoms
            #-------------------------------------
            # Keep only small distances
            distanceatm=cell_mod.distance( atoms1, cell2, atom1, atom2 )
            if distanceatm < 3.0
                #_-------------------------------
                # Move on of the atom in the closest image to the other
                #--------------------------------------
                # Keep orignal position in mind
                temp1=atoms1.positions[atom1,:]
                temp2=atoms1.positions[atom2,:]
                atoms.positions[atom1,:]-=ELF1.origin
                atoms.positions[atom2,:]-=ELF1.origin
                # Moving atom2
                for j=1:3
                    di = atoms1.positions[atom1,j]-atoms1.positions[atom2,j]
                    if di > cell2.length[j]*0.5
                        atoms1.positions[atom2,j] += cell2.length[j]
                    end
                    if di < -cell2.length[j]*0.5
                        atoms1.positions[atom2,j] -= cell2.length[j]
                    end
                end
                #--------------------------------------
                # Now we have the position of the center...
                position_origin=(atoms1.positions[atom1,:]+atoms1.positions[atom2,:])/2
                # Wrapping it in PBC
                for j=1:3
                    position_origin[j] = cell_mod.wrap( position_origin[j] , cell2.length[j] )
                end
                #-----------------------------------------
                # First guess of index
                #--------------------------------------
                floats=[0,0,0]
                for j=1:3
                    check=position_origin[j]*ELF1.nb_vox[j]/cell2.length[j]
                    floats[j]=trunc(check)
                    if check - floats[j] > 0.5
                        floats[j] += 1
                    end
                    if floats[j] > ELF1.nb_vox[j]-1
                        floats[j]=0
                    end
                end
                #-----------------------------------------
                # Correct the position...
                distance1=0
                for l=1:3
                    distance1+=cell_mod.dist1D( floats[l]*cell2.length[l]/ELF1.nb_vox[l], position_origin[l], cell2.length[l] )^2
                end
                for p=-1:1:1
                    for j=-1:1:1
                        for k=-1:1:1
                            if k != 0 || j != 0 || p != 0
                                new_floats=[0,0,0]
                                new_floats[1]=floats[1]+p
                                new_floats[2]=floats[2]+j
                                new_floats[3]=floats[3]+k
                                for l=1:3
                    		    		if new_floats[l] < 0
                    					new_floats[l] = ELF1.nb_vox[l] - new_floats[l]
                    		    		end
                    		    		if new_floats[l] >= ELF1.nb_vox[l]-1
                    					new_floats[l] = new_floats[l]-(ELF1.nb_vox[l]-1)
                    		    		end
                    			end
                                distance2=0
                                for l=1:3
                                    distance2+=cell_mod.dist1D( new_floats[l]*cell2.length[l]/ELF1.nb_vox[l], position_origin[l], cell2.length[l] )^2
                                end
                                if distance2 < distance1
                                    distance1 = distance2
                                    floats=new_floats
                                end
                            end
                        end
                    end
                end
                for i=1:3
                    floats[i] += 1
                end
                write(file_check,string(distanceatm," ",distance1," ",ELF1.matrix[ floats[1], floats[2], floats[3]]," ",sqrt(3)/2*cell2.length[1]/ELF1.nb_vox[1], "\n"))
                elf_value=ELF1.matrix[ floats[1], floats[2], floats[3]]
                if elf_value < 0.7 && distanceatm < 1.1
                    write(log,string(i," ",atom1," ",atom2," ",elf_value," ",distanceatm,"\n"))
                end
                #------------------------------------------------------------------------------------------------------------------------

                # Histogram
                #------------------------------------------------------------------------------------
                if elf_value  < max_elf && distanceatm < max_dist
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
                #------------------------------------------------------------------------------------
                # Putting the position back to normal for futher analysis
                atoms1.positions[atom1,:]=temp1
                atoms1.positions[atom2,:]=temp2
            end
        end
    end
end

close(file_check)
close(log)

filehist=open("/home/moogmt/hist_elf_GLYCINE_shoot1_OH.dat","w")
for i=1:nelf
    for j=1:ndist
        write(filehist,string( i*ddist+min_dist, " ", j*delf+min_elf," ",hist2d[i,j]/count2,"\n"))
    end
    write(filehist,"\n")
end
close(filehist)
