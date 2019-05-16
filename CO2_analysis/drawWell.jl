GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))


folder_base="/home/moogmt/"

N_points=500
max_=1
min_=0
I0=[9,15,7]
x0=[0.15,0.5,0.8]
del=[0.1,0.2,0.05]
delta=(max_-min_)/N_points
file_out=open(string(folder_base,"wells.dat"),"w")
for i=1:N_points
    value=0
    for j=1:3
        value += -I0[j]*exp(-(i*delta-x0[j])*(i*delta-x0[j])/(2*del[j]*del[j]))
    end
    write(file_out,string(i*delta," ",value,"\n"))
end
close(file_out)
