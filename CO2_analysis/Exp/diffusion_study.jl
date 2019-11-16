GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using exp_data

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

func="PBE-MT"

unit=0.005
nbC = 32
nbO = 64

V=8.82
T=3000

folder=string(folder_base,V,"/",T,"K/")
folder_out=string(folder,"Data/")
file="TRAJEC.xyz"

traj = filexyz.readFastFile( string(folder,file) )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size( traj )[1]
nb_atoms=size( traj[1].names )[1]


barycenter=zeros(Real,nb_steps,3)
C_barycenter=zeros(Real,nb_steps,3)
O_barycenter=zeros(Real,nb_steps,3)
d_barycenter=zeros(Real,nb_steps)
d_barycenterC=zeros(Real,nb_steps)
d_barycenterO=zeros(Real,nb_steps)

file_out=open(string(folder_out,"barycenter_diff.dat"),"w")
for step=1:nb_steps
    print("Progress:",step/nb_steps*100,"\n")
    for i=1:3
        for carbon=1:nbC
            barycenter[step,i] += traj[step].positions[carbon,i]
            C_barycenter[step,i] += traj[step].positions[carbon,i]
        end
        C_barycenter[step,i] /= nbC
        for oxygen=1:nbO
            barycenter[step,i] += traj[step].positions[nbC+oxygen,i]
            O_barycenter[step,i] += traj[step].positions[nbC+oxygen,i]
        end
        O_barycenter[step,i] /= nbO
        barycenter[step,i] /= nb_atoms
    end
    for i=1:3
        d_barycenter[step] += (barycenter[step,i]-barycenter[1,i])*(barycenter[step,i]-barycenter[1,i])
        d_barycenterC[step] += (C_barycenter[step,i]-C_barycenter[1,i])*(C_barycenter[step,i]-C_barycenter[1,i])
        d_barycenterO[step] += (O_barycenter[step,i]-O_barycenter[1,i])*(O_barycenter[step,i]-O_barycenter[1,i])
    end
    d_barycenter[step] = sqrt(d_barycenter[step])
    d_barycenterC[step] = sqrt(d_barycenterC[step])
    d_barycenterO[step] = sqrt(d_barycenterO[step])
    write(file_out,string(step," ",d_barycenter[step]," ",d_barycenterC[step]," ",d_barycenterO[step],"\n"))
end
close(file_out)

barycenter=[]
C_barycenter=[]
O_barycenter=[]

file_position_C=open(string(folder_out,"positions_C.dat"),"w")
for step=1:nb_steps
    write(file_position_C,string(step," "))
    for carbon=1:nbC
        for i=1:3
            write(file_position_C,string(traj[step].positions[carbon,i]," "))
        end
    end
    write(file_position_C,string("\n"))
end
close(file_position_C)

file_position_O=open(string(folder_out,"positions_O.dat"),"w")
for step=1:nb_steps
    write(file_position_O,string(step," "))
    for oxygen=nbC+1:nb_atoms
        for i=1:3
            write(file_position_O,string(traj[step].positions[oxygen,i]," "))
        end
    end
    write(file_position_O,string("\n"))
end
close(file_position_O)

min_x=0
max_x=1.5
nb_box=100
delta=(max_x-min_x)/nb_box
hist1D=zeros(Real,100)
for carbon=1:nbC
    print("Progress: ",carbon/nbC*100,"%\n")
    for step=2:nb_steps
        for j=1:nb_box
            distance=0
            for k=1:3
                distance+=(traj[step-1].positions[carbon,k]-traj[step].positions[carbon,k])*(traj[step-1].positions[carbon,k]-traj[step].positions[carbon,k])
            end
            distance=sqrt(distance)
            if  distance > (j-1)*delta && distance < j*delta
                hist1D[j] += 1
                break
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"histDX-carbon.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta," ",hist1D[i],"\n"))
end
close(file_out)

cut_off_distance=0.3
nb_points=nb_steps
rho=zeros(Real,nb_points,nbC)
for carbon=1:nbC
    print("Progress: ",carbon/nbC*100,"%\n")
    for i=1:nb_points
        for j=i-500:i+500
            if j < 1 || j > nb_steps
                continue
            end
            distance=0
            for k=1:3
                distance += (traj[i].positions[carbon,k]-traj[j].positions[carbon,k])*(traj[i].positions[carbon,k]-traj[j].positions[carbon,k])
            end
            if distance < cut_off_distance
                rho[i,carbon] += 1
            end
        end
    end
end

min_x=0
max_x=1000
nb_box=200
delta=(max_x-min_x)/nb_box
hist1D=zeros(Real,nb_box)
for carbon=1:nbC
    print("Progress: ",carbon/nbC*100,"%\n")
    for i=1:nb_points
        for j=1:nb_box
            if rho[i,carbon] > (j-1)*delta && rho[i,carbon] < j*delta
                hist1D[j] += 1
            end
        end
    end
end
hist1D/=sum(hist1D)

file_out=open(string(folder_out,"hist-carbon.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*delta," ",hist1D[i],"\n"))
end
close(file_out)
