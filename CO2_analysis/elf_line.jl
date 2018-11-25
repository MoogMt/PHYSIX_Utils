include("atoms.jl")
include("cell.jl")
include("cubefile.jl")


step_min=1
step_max=1001
d_step=10

start_distance=1.2
stoping_distance=2.0
nb_points=25

distance_data=[]
elf_data=[]
used=[]

for step=step_min:d_step:step_max

    #------------------------------------------------------------------------------
    atoms, cell1, ELF1 = cube_mod.readCube(string("/media/moogmt/Stock/cube_TS14_gly/ELF_shoot1_",step,".cube"))
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

    for atom1=3:83
        for atom2=85:size(atoms.names)[1]
            #------------------------------------------------------------------------------
            distanceatm=cell_mod.distance( atoms, cell2, atom1, atom2 )
            #------------------------------------------------------------------------------

            if distanceatm > start_distance && distanceatm < stoping_distance

                #-----------------
                # Ajusting atom2
                #------------------------------------------------------------------------------
                temp21=atoms.positions[atom2,:]
                for j=1:3
                    di = atoms.positions[atom1,j]-temp21[j]
                    if di > cell2.length[j]*0.5
                        temp21 += cell2.length[j]
                    end
                    if di < -cell2.length[j]*0.5
                        temp21 -= cell2.length[j]
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
                v=temp21 - atoms.positions[atom1,:]
                #------------------------------------------------------------------------------

                #------------------------------------------------------------------------------
                points=zeros(nb_points,3)
                for i=1:nb_points
                    points[i,:]=atoms.positions[atom1,:]+v/nb_points*i
                end
                #------------------------------------------------------------------------------

                #--------------
                # Wrapping PBC
                #------------------------------------------------------------------------------
                for j=1:nb_points
                    for i=1:3
                        points[j,i]=cell_mod.wrap(points[j,i], params[i])
                    end
                end
                #------------------------------------------------------------------------------

                #------------------------------------------------------------------------------
                index=zeros(nb_points,3)
                for j=1:nb_points
                    for i=1:3
                        check=points[j,i]*ELF1.nb_vox[i]/params[i]
                        index[j,i]=trunc(check)
                        if check - index[j,i] > 0.5
                            index[j,i] += 1
                        end
                        if index[j,i] > ELF1.nb_vox[i]-1
                            index[j,i]=0
                        end
                    end
                end
                #------------------------------------------------------------------------------

                #------------------------------------------------------------------------------
                distance1=zeros(nb_points)
                for j=1:nb_points
                    for i=1:3
                        distance1[j]+=cell_mod.dist1D( index[j,i]*cell2.length[i]/ELF1.nb_vox[i], points[j,i], cell2.length[i] )^2
                    end
                end
                #------------------------------------------------------------------------------

                #------------------------------------------------------------------------------
                for j=1:nb_points
                    for l=-1:1:1
                        for m=-1:1:1
                            for n=-1:1:1
                                new_index=[0,0,0]
                                new_index[1]=index[j,1]+l
                                new_index[2]=index[j,2]+m
                                new_index[3]=index[j,3]+n
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
                                    distance2+=cell_mod.dist1D( new_index[i]*cell2.length[i]/ELF1.nb_vox[i], points[j,i], params[i] )^2
                                end
                                if distance2 < distance1[j]
                                    distance1[j] = distance2
                                    index[j,:]=new_index
                                end
                            end
                        end
                    end
                end
                #------------------------------------------------------------------------------

                #------------------
                index+=1
                #------------------

                #-------------
                # Adding data
                #---------------------------------------------------
                for j=1:nb_points
                    push!( used, 0 )
                    push!( distance_data, j )
                    push!( elf_data, ELF1.matrix[ Int(index[j,1]), Int(index[j,2]), Int(index[j,3]) ] )
                end
                #---------------------------------------------------

                #endif
            end
            # end atom2
        end
        # end atom1
    end
    print("step:",step,"\n")
    # end step
end

min_elf=0
max_elf=1.0
n_elf=100
d_elf=(max_elf-min_elf)/n_elf

hist_2D=zeros(n_elf,nb_points)

count=0

for i=1:n_elf
    print("progress: ",i/n_elf*100,"%\n")
    for j=1:nb_points
        for k=1:size(used)[1]
            if used[k] == 0
                if i*d_elf+min_elf < elf_data[k] && (i+1)*d_elf+min_elf > elf_data[k] && j == distance_data[k]
                    used[k] += 1
                    hist_2D[i,j] += 1
                    count += 1
                end
            end
        end
    end
end

hist_2D /= count

file_hist=open(string("/home/moogmt/lineOH_elf_bond_start",start_distance,"_stop_",stoping_distance,".dat"),"w")
for i=1:n_elf
    for j=1:nb_points
        write(file_hist,string(j," ",i*d_elf+min_elf," ",hist_2D[i,j],"\n"))
    end
    write(file_hist,"\n")
end
close(file_hist)
