GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion

# Folder for data


folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2


T=3000
V=8.82

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

# if ! isfile(file)
#     continue
# end

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nbC=32
nbO=64

nb_steps=20000

rmin=0
rmax=V/2
dr=0.01
nb_box=Int((rmax-rmin)/dr)
gr=zeros(nb_box)

for step=1:nb_steps
    for carbon=1:nbC
        for oxygen=1:nbO
            dist=cell_mod.distance(traj[step],cell,carbon,nbC+oxygen)
            if dist < rmax
                gr[ Int( trunc( (dist-rmin)/dr ) )  ] += 1
            end
        end
    end
end


file_out=open(string(folder_out,"gr_real.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*dr+rmin," ",gr[i]/(4*pi*(i*dr+rmin)^2),"\n"))
end
close(file_out)

rho=2

n_sk=1000
sk=ones(n_sk)
dk=0.05

for i=1:n_sk
    int_=0
    for j=1:nb_box
        int_ += (j*dr+rmin)^2*gr[j]*sin((i*dk*(j*dr+rmin)))/((i*dk*(j*dr+rmin)))*dr
    end
    sk[i] += 4*pi*rho*int_
end

file_out=open(string(folder_out,"sk_real.dat"),"w")
for i=1:nb_box
    write(file_out,string(i*dk," ",sk[i],"\n"))
end
close(file_out)

y=zeros(20000)
x=zeros(20000)
dx=0.05

for i=1:1000
    x[i] = i*dx
    y[i] =sin(i*dx)
end
ddx=dct(x)
ddy=dct(y)

file_out=open(string(folder_out,"test_truc.dat"),"w")
for i=1:nb_box
    write(file_out,string(x[i]," ",y[i]," ",ddx[i]," ",ddy[i]," \n"))
end
close(file_out)
