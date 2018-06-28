include("contactmatrix.jl")
using PyPlot


volumes=[9.8,9.4,9.3,9.2,9.1,9.0,8.82]

T=3000
for V in volumes

folder=string("/home/moogmt/CO2/CO2_AIMD/",V,"/",T,"K/")
file=string("gr_ALL_",V,"_",T,"K.dat")

file2=open(string(folder,file))
lines=readlines(file2)
close(file2)

nb_points=size(lines)[1]

x=zeros(nb_points-1)
g_r=zeros(nb_points-1)

for i=1:nb_points-1
    line=split(lines[i])
    x[i]=parse(Float64,line[1])
    g_r[i]=parse(Float64,line[2])
end
rho=1.

s_q =abs(fft(g_r))
q=abs(x)

file3=open(string("/home/moogmt/CO2/CO2_AIMD/sq_",V,"_",T,"K.dat"),"w")
for i=1:nb_points-1
    write(file3,string(q[i]," ",s_q[i],"\n"))
end
close(file3)

end
