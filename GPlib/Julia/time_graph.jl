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

function average(x)
    avg=0
    for i=1:size(x)[1]
        avg += x[i]
    end
    return avg/size(x)[1]
end

# Current Volume and Temperature
current_volume=parse(Float64,"9.4")
# for T in [2000,2500,3000]
T=3000
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",current_volume,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
unit=0.0005
stride=1

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_atoms=size(atoms[1].names)[1]

cutoff=1.7


nb_steps=size(atoms)[1]
steps=nb_steps/10
time_corr_tot=zeros(steps)

for atom=1:32

bond_matrix=zeros(size(atoms)[1],64)
for step1=1:size(atoms)[1]
    for i=33:96
        if i != atom
            if cell_mod.distance(atoms[step1],cell,i,atom) < 1.7
                bond_matrix[step1,i-32] += 1
            end
        end
    end
end

print("progress:",atom/32*100,"%\n")

time_corr=zeros(steps)
for i=33:nb_atoms
    time_corr_loc=autocorrelation(bond_matrix[:,i-32],0,steps,steps,nb_steps-steps)
    for j=1:steps
        time_corr[j] += time_corr_loc[j]
    end
end

for i=1:steps
    time_corr_tot[i] += time_corr[i]
end

end

value=time_corr_tot[1]
for i=1:steps
    time_corr_tot[i] = time_corr_tot[i]/value
end

file_check=open(string("/home/moogmt/time_",current_volume,"-",T,".dat"),"w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",time_corr_tot[Int(i)],"\n"))
end
close(file_check)
# end
