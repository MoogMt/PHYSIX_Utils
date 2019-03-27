GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"cell.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5


# cols (both files):
#   0:   natoms x nfi (natoms x 1, natoms x 2, ...)
#   1-3: x,y,z cartesian coords [Bohr]
#   4-6: x,y,z cartesian velocites [Bohr / thart ]
#        thart = Hartree time =  0.024189 fs
# FTRAJECTORY extra:
#   7-9: x,y,z cartesian forces [Ha / Bohr]

V=8.82
T=3000


file_in=open(string(folder_base,V,"/",T,"K/FTRAJECTORY2"))
folder_out=string(folder_base,V,"/",T,"K/")
lines=readlines(file_in)
close(file_in)

nb_atoms=96
nb_steps=500

cell=cell_mod.Cell_param(V,V,V)

positions=zeros(Real,nb_steps,nb_atoms,3)
forces=zeros(Real,nb_steps,nb_atoms,3)

start_step=2500

for i=start_step+1:start_step+nb_steps
    for j=1:nb_atoms
        line=split(lines[(i-1)*nb_atoms+j])
        for k=1:3
            positions[i-start_step,j,k] = cell_mod.wrap(parse( Float64, line[k+1] ),cell.length[k])
            forces[i-start_step,j,k]    = parse( Float64, line[k+7] )
        end
    end
end

positions_points=zeros(Real,nb_max_neigh,3)
for atom1=1:nb_max_neigh
    for atom2=1:nb_atoms
    distance=cell_mod.distance(positions[1,1,:],positions[1,nbC+oxygen,:],cell.length[:])
    end
end

theta_value=100
phi_value=200
r=1.5

matrix=zeros(theta_value,phi_value)*
point=zeros(Real,3)
for i=1:theta_value
    theta=i*pi/theta_value
    for j=1:phi_value
        phi=j*2*pi/phi_value
        point[1]=r*cos(theta)*cos(phi)
        point[2]=r*cos(theta)*sin(phi)
        point[3]=r*sin(theta)
        for j=1:nb_max_neigh

        end
    end
end


f_min=0
f_max=0
for step=1:nb_steps
    for atom=1:nb_atoms
        for k=1:3
            if forces[step,atom,k] > f_max
                global f_max = forces[step,atom,k]
            end
            if forces[step,atom,k] < f_min
                global f_min = forces[step,atom,k]
            end
        end
    end
end

nb_box=1000
delta=(f_max-f_min)/nb_box

hist1D=zeros(nb_box)
for step=1:nb_steps
    for atom=1:nb_atoms
        for k=1:3
            for i=1:nb_box
                if (i-1)*delta+f_min < forces[step,atom,k] && i*delta+f_min > forces[step,atom,k]
                    hist1D[i] += 1
                    break
                end
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"hist_forces_comp.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta+f_min," ",hist1D[i],"\n"))
end
close(file_out)

f_min=100
f_max=0
forces_norm=zeros(nb_steps,nb_atoms)
for step=1:nb_steps
    for atom=1:nb_atoms
        for k=1:3
            forces_norm[step,atom] += forces[step,atom,k]*forces[step,atom,k]
        end
        forces_norm[step,atom] = sqrt(forces_norm[step,atom])
        if forces_norm[step,atom] > f_max
            global f_max=forces_norm[step,atom]
        end
        if forces_norm[step,atom] < f_min
            global f_min=forces_norm[step,atom]
        end
    end
end

nb_box=1000
delta=(f_max-f_min)/nb_box

hist1D=zeros(nb_box)
for step=1:nb_steps
    for atom=1:nb_atoms
        for k=1:3
            for i=1:nb_box
                if (i-1)*delta+f_min < forces[step,atom,k] && i*delta+f_min > forces[step,atom,k]
                    hist1D[i] += 1
                    break
                end
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"hist_forces_norm.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta+f_min," ",hist1D[i],"\n"))
end
close(file_out)


f_minC=100
f_maxC=0
f_minO=100
f_maxO=0
forces_C=zeros(nb_steps,nbC)
forces_O=zeros(nb_steps,nbO)
for step=1:nb_steps
    for atom=1:nbC
        for k=1:3
            forces_C[step,atom] += forces[step,atom,k]*forces[step,atom,k]
        end
        forces_C[step,atom] = sqrt(forces_C[step,atom])
        if forces_C[step,atom] > f_maxC
            global f_maxC=forces_C[step,atom]
        end
        if forces_C[step,atom] < f_minC
            global f_minC=forces_C[step,atom]
        end
    end
    for k=1:3
    for atom=1:nbO
            forces_O[step,atom] += forces[step,nbC+atom,k]*forces[step,nbC+atom,k]
        end
        forces_O[step,atom] = sqrt(forces_O[step,atom])
        if forces_O[step,atom] > f_maxO
            global f_maxO=forces_O[step,atom]
        end
        if forces_O[step,atom] < f_minO
            global f_minO=forces_O[step,atom]
        end
    end
end

nb_box=1000
delta=(f_maxC-f_minC)/nb_box

hist1D=zeros(nb_box)
for step=1:nb_steps
    for atom=1:nbC
        for i=1:nb_box
            if (i-1)*delta+f_min < forces_C[step,atom] && i*delta+f_min > forces_C[step,atom]
                hist1D[i] += 1
                break
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"hist_forcesC_norm.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta+f_min," ",hist1D[i],"\n"))
end
close(file_out)

nb_box=1000
delta=(f_maxO-f_minO)/nb_box

hist1D=zeros(nb_box)
for step=1:nb_steps
    for atom=1:nbC
        for i=1:nb_box
            if (i-1)*delta+f_min < forces_O[step,atom] && i*delta+f_min > forces_O[step,atom]
                hist1D[i] += 1
                break
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"hist_forcesO_norm.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta+f_min," ",hist1D[i],"\n"))
end
close(file_out)





max_neigh=10
distancesO=zeros(nbC*nb_steps,max_neigh)
forcesO=zeros(nbC*nb_steps,max_neigh)
for i=1:nb_steps
    for j=1:nbC
        local_distance=zeros(nbC)
        forces_local=zeros(nbC)
        for k=1:nbO
            for l=1:3
                local_distance[l] += (positions[i,j,l]-positions[i,nbC+k,l])^2
            end
            local_distance[l] = sqrt( local_distance[l] )
            forces_local[i,j,k]
        end
    end
end
