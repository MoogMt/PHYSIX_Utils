include("atoms.jl")
include("cell.jl")
include("cubefile.jl")

atom1=20
atom2=62

file=open(string("/home/moogmt/check",atom1,"_",atom2,".dat"),"w")

nb_pt_line=30

current_volume=8.82
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

file_test=open(string("/home/moogmt/check",atom1,"_",atom2,".xyz"),"w")
# write(file_test,string("150\nCHECK\n"))
file_atom=open("/home/moogmt/test_line.xyz","w+")
for i=0:49
    print(string(i,"\n"))
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82_dyn/",i,"_structure/")
    atoms1, cell1, ELF1 = cube_mod.readCube(string(folder,"ELF.cube"))
    # PBC for orthorombic structure
    for j=1:3
        di = atoms1.positions[atom1,j]-atoms1.positions[atom2,j]
        if di > cell.length[j]*0.5
            atoms1.positions[atom1,j] -= cell.length[j]
        end
        if di < cell.length[j]*0.5
            atoms1.positions[atom1,j] += cell.length[j]
        end
    end
    pos1=cube_mod.getClosest(( atoms1.positions[atom1,:]+atoms1.positions[atom2,:])/2 ,ELF1)
    write(file,string( cell_mod.distance(atoms1,cell,atom1,atom2), " ", ELF1.matrix[ pos1[1], pos1[2], pos1[3] ] ," ",i, "\n") )
    # file_temp=open(string("/home/moogmt/line_",i,"_",atom1,"-",atom2,".dat"),"w")
    # Vector to trace line
    #========================#
    dv= (atoms1.positions[atom2,:]-atoms1.positions[atom1,:])/nb_pt_line
    norm=0
    for j=1:3
        norm += dv[j]^2
    end
    norm=sqrt(norm)
    #========================#
    if i== 0
        write(file_atom,string(nb_pt_line+2,"\nCHECK\n"))
        write(file_atom,"C ")
        for j=1:3
            write(file_atom,string(atoms1.positions[atom1,j]," "))
        end
        write(file_atom,"\n")
        write(file_atom,string("O "))
        for j=1:3
            write(file_atom,string(atoms1.positions[atom2,j]," "))
        end
        write(file_atom,"\n")
    end
    #==========================#
    for j=1:nb_pt_line
        posline=cube_mod.getClosest( atoms1.positions[atom1,:] + j*dv ,ELF1)
        write(file_temp,string(j*norm," ",ELF1.matrix[ posline[1], posline[2], posline[3] ],"\n"))
        if i==0
        write(file_atom,string("N "))
        for k=1:3
            write(file_atom,string(posline[k]*8.82/60," "))
        end
        write(file_atom,"\n")
        end
    end
    #==========================#
    # close(file_temp)
    # write(file_test,string("C "))
    # for j=1:3
    #     write(file_test,string(atoms1.positions[atom1,j]," "))
    # end
    # write(file_test,"\n")
    # write(file_test,string("O "))
    # for j=1:3
    #     write(file_test,string(atoms1.positions[atom2,j]," "))
    # end
    # write(file_test,"\n")
    # write(file_test,string("N "))
    # for j=1:3
    #     write(file_test,string(pos1[j]*8.82/60," "))
    # end
    # write(file_test,"\n")
end
close(file_test)

close(file)
close(file_atom)
