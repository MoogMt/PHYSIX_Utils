folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/9.8/3000K/"

file="FTRAJECTORY"

nb_step=0

file=string(folder,file)


print(nb_step)

file_count=open(file)
while !eof(file_count)
    print("check")
end
close(file_count)

velocities=zeros(nb_steps,3)
vnorms=zeros(nb_steps,1)

file_traj=open(file)
for i=1:nb_steps
    line=split(readline(file_traj))
    for j=1:3
        velocities[i,j]=parse(Float64,line[j])
        vnorms[j] += velocities[i,j]^2
    end
    vnorms[j] = sqrt(vnorms[j])
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

data=fft(autocorrelation2(vnorms,frac))

using PyPlot

plot(data,"r-.")
