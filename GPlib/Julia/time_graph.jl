include("contactmatrix.jl")

function autocorrelation(x,avg,N,i,M)
    C=zeros(N)
    for k=i:i+M
        for n=1:N-1
            C[n] += (x[k]-avg)*(x[k-n]-avg)
        end
    end
    return C
end

function autocorrelation2()
end

function average(x)
    avg=0
    for i=1:size(x)[1]
        avg += x[i]
    end
    return avg/size(x)[1]
end

# Current Volume and Temperature

Volumes=[9.4]

for volume in Volumes
# for T in [2000,2500,3000]
T=3000
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",volume,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
unit=0.0005*5
stride=1
frac=0.8
cutoff=1.7

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(volume,volume,volume)

nb_atoms=size(atoms[1].names)[1]
nb_steps=size(atoms)[1]
steps=Int(trunc(nb_steps*frac))

time_corr=zeros(steps)

for atom1=1:32
    print("progress:",atom1/32*100,"%\n")
    for atom2=33:nb_atoms

        bond_matrix=zeros(nb_steps)
        for step=1:nb_steps
            if cell_mod.distance(atoms[step],cell,atom1,atom2) < 1.7
                bond_matrix[step] = 1
            end
        end

        time_corr +=autocorrelation(bond_matrix,0,steps,steps,nb_steps-steps)

    end
end

# Normalization
value=time_corr[1]
for i=1:steps
    time_corr[i] = time_corr[i]/value
end

file_check=open(string("/home/moogmt/time_",volume,"-",T,".dat"),"w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",time_corr[Int(i)],"\n"))
end
close(file_check)
# end
end
