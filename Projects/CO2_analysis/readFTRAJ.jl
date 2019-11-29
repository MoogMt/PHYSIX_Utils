GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Thermo data
# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5



V=8.82
T=3000

folder_out=string("folder_base",V,"/",T,"K/")

file_in=open(string(folder_base,V,"/",T,"K/FTRAJECTORY2"))
lines=readlines(file_in)
close(file_in)

nb_atoms=96
nb_steps=Int(size(lines)[1]/96)

cell=cell_mod.Cell_param(V,V,V)

positions=zeros(Real,nb_steps,nb_atoms,3)
forces=zeros(Real,nb_steps,nb_atoms,3)

start_step=0

# parsing FTRAJECTORY file
for i=start_step+1:start_step+nb_steps
    for j=1:nb_atoms
        line=split(lines[(i-1)*nb_atoms+j])
        for k=1:3
            positions[i-start_step,j,k] = cell_mod.wrap(parse( Float64, line[k+1] ),cell.length[k])
            forces[i-start_step,j,k]    = parse( Float64, line[k+7] )
        end
    end
end

#----------------------------------------#

max_force=0
min_force=5000
for step=1:nb_steps
    for atom=1:nb_atoms
        force_norm=0
        for i=1:3
            force_norm += forces[step,atom,i]*forces[step,atom,i]
        end
        force_norm = sqrt(force_norm)
        if force_norm > max_force
            global max_force = force_norm
        end
        if force_norm < min_force
            global min_force = force_norm
        end
    end
end

nb_box=200
delta_box=(max_force-min_force)/nb_box
hist_force_norm=zeros(Real,nb_box+1)
for step=1:nb_steps
    for atom=1:nb_atoms
        force_norm=0
        for i=1:3
            force_norm += forces[step,atom,i]*forces[step,atom,i]
        end
        hist_force_norm[ Int(trunc((sqrt(force_norm)-min_force)/delta_box + 1)) ] += 1
    end
end
hist_force_norm /= sum(hist_force_norm)

file_out=open(string(folder_out,"Force_Hist.dat"),"w")
for i=1:nb_box+1
    write(file_out,string(i*delta_box+min_force," ",hist_force_norm[i],"\n"))
end
close(file_out)
#-----------------------------------------------#

max_force=0
min_force=5000
for step=1:nb_steps
    for atom=1:nbC
        force_norm=0
        for i=1:3
            force_norm += forces[step,atom,i]*forces[step,atom,i]
        end
        force_norm = sqrt(force_norm)
        if force_norm > max_force
            global max_force = force_norm
        end
        if force_norm < min_force
            global min_force = force_norm
        end
    end
end

nb_box=200
delta_box=(max_force-min_force)/nb_box
hist_force_norm=zeros(Real,nb_box+1)
for step=1:nb_steps
    for atom=1:nbC
        force_norm=0
        for i=1:3
            force_norm += forces[step,atom,i]*forces[step,atom,i]
        end
        hist_force_norm[ Int(trunc((sqrt(force_norm)-min_force)/delta_box + 1)) ] += 1
    end
end
hist_force_norm /= sum(hist_force_norm)

file_out=open(string(folder_out,"Force_Hist_C.dat"),"w")
for i=1:nb_box+1
    write(file_out,string(i*delta_box+min_force," ",hist_force_norm[i],"\n"))
end
close(file_out)

#-----------------------------------------------#

max_force=0
min_force=5000
for step=1:nb_steps
    for atom=1:nbO
        force_norm=0
        for i=1:3
            force_norm += forces[step,nbC+atom,i]*forces[step,nbC+atom,i]
        end
        force_norm = sqrt(force_norm)
        if force_norm > max_force
            global max_force = force_norm
        end
        if force_norm < min_force
            global min_force = force_norm
        end
    end
end

nb_box=200
delta_box=(max_force-min_force)/nb_box
hist_force_norm=zeros(Real,nb_box+1)
for step=1:nb_steps
    for atom=1:nbO
        force_norm=0
        for i=1:3
            force_norm += forces[step,nbC+atom,i]*forces[step,nbC+atom,i]
        end
        hist_force_norm[ Int(trunc((sqrt(force_norm)-min_force)/delta_box + 1)) ] += 1
    end
end
hist_force_norm /= sum(hist_force_norm)

file_out=open(string(folder_out,"Force_Hist_O.dat"),"w")
for i=1:nb_box+1
    write(file_out,string(i*delta_box+min_force," ",hist_force_norm[i],"\n"))
end
close(file_out)
