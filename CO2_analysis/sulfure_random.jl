GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")


include(string(CO2folder,"markovCO2.jl"))


# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]

# Cut-off distance for bonds
cut_off_bond = 1.75

# Volume of the cell (only orthorombic is implemented yet)


# Defining base folder
folder_in=string("/media/moogmt/Stock/Mathieu/Sulfur/")

V=16.0

if ! ispath( string(folder_in,"/",V,"/") )
    mkdir(string(folder_in,"/",V,"/"))
end

# Folder where to put the out data
folder_out=string(folder_in,"/",V,"/")

nb_atoms=96
cut_off = 2.0

for random_conf=1:20

positions_confirmed=zeros(nb_atoms,3)

count_atom=0
while count_atom < 96
    pos=rand(3)*V
    accepted=true
    for atom=1:count_atom
        dist=0
        for i=1:3
            dist+= cell_mod.dist1D(positions_confirmed[atom,i],pos[i],V)^2
        end
        if sqrt(dist) < cut_off
            accepted=false
            break
        end
    end
    if accepted
        count_atom += 1
        positions_confirmed[count_atom,:] = pos[:]
    end
end

file_out=open(string(folder_out,"S-",random_conf,".cpmd"),"w")
for i=1:nb_atoms
    for j=1:3
        write(file_out,string(positions_confirmed[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)

file_out=open(string(folder_out,"S-",random_conf,".xyz"),"w")
write(file_out,string(nb_atoms,"\n"))
write(file_out,string("TEST \n"))
for i=1:nb_atoms
    write(file_out,string("S "))
    for j=1:3
        write(file_out,string(positions_confirmed[i,j]," "))
    end
    write(file_out,string("\n"))
end
close(file_out)
end
