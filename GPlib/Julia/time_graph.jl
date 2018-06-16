include("contactmatrix.jl")

function autocorrelation(x,N,i,M)
    C=zeros(N)
    for k=i:i+M-1
        for n=1:N-1
            C[n] += x[k]*x[k-n]
        end
    end
    return C
end

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

function average(x)
    avg=0
    for i=1:size(x)[1]
        avg += x[i]
    end
    return avg/size(x)[1]
end

# Current Volume and Temperature

Volumes=[8.82]

#for volume in Volumes
# for T in [2000,2500,3000]
#T=3000
folder=string("/home/moogmt/8.82/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
unit=0.0005*5
stride=4
frac=0.8
cutoff=1.7

volume=8.82

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(volume,volume,volume)

nb_atoms=size(atoms[1].names)[1]
nb_steps=size(atoms)[1]
steps=Int(trunc(nb_steps*frac))

time_corr=zeros(steps)

max=0
corr_stock=zeros(steps)
index1=0
index2=0

for atom1=1:32
    print("progress:",atom1/32*100,"%\n")
    for atom2=33:nb_atoms
        bond_matrix=zeros(nb_steps)
        for step=1:nb_steps
            if cell_mod.distance(atoms[step],cell,atom1,atom2) < 1.7
                bond_matrix[step] = 1
            end
        end
        corr=autocorrelation2(bond_matrix,frac)
        sum=0
        for i=1:size(corr)[1]
            sum+=corr[i]
        end
        if sum > max
            corr_stock = corr
            max = sum
            index1=atom1
            index2=atom2
        end
        time_corr += corr
    end
end

# Normalization
value=time_corr[1]
time_corr /= value

value=corr_stock[1]
corr_stock /= value

print("atom1:",index1," atom2: ",index2,"\n")

file_check=open(string("/home/moogmt/time_long.dat"),"w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",corr_stock[Int(i)],"\n"))
end
close(file_check)

file_check=open(string("/home/moogmt/time.dat"),"w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",time_corr[Int(i)],"\n"))
end
close(file_check)
# end
#end

file_dist=open(string("/home/moogmt/dist.dat"),"w")
for step=1:nb_steps
    write(file_dist,string(step*stride*unit," ",cell_mod.distance(atoms[step],cell,8,42),"\n"))
end
close(file_dist)
