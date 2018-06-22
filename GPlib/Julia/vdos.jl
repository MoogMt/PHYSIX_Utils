folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/9.8/3000K/"

file="FTRAJECTORY"

file2=string(folder,file)
print(nb_step)

step = 0
step +=1
open( file2 ) do f
    while !eof(f)
        l=readline(f)
        step = step + 1
    end
end

nb_steps=10000
velocities=zeros(nb_steps,96,3)
vnorms=zeros(nb_steps,96)

file_traj=open(file2)
for i=1:nb_steps
    line=split(readline(file_traj))
    for k=1:96
        for j=1:3
            velocities[i,k,j]=parse(Float64,line[j])
            vnorms[j,k] += velocities[i,k,j]^2
        end
        vnorms[i,k] = sqrt(vnorms[i,k])
    end
end
close(file_traj)

function autocorrelation2(x,frac)
    M=size(x)[1]
    N=Int(trunc(M*frac))
    C=zeros(N)
    for n=1:N
        for m=1:M-n
            C[n] += x[m]*x[m+n]
        end
        C[n] = C[n]/(M-n)
    end
    return C
end

frac=0.3

data=zeros(nb_steps*frac)

for i=1:96
    data[:] += autocorrelation2(vnorms[:,i],frac)
end
data /= 96

dataplot=fft(data)

using PyPlot

plot(data,"r-.")
