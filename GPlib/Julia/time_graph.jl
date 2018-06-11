include("contactmatrix.jl")

# Current Volume and Temperature
current_volume=parse(Float64,"8.82")
for T in [2000,2500,3000]
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",current_volume,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
unit=0.0005*5
stride=1

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

cutoff=1.7
atom=1
bond_matrix=zeros(nb_steps,Int(trunc((nb_atoms*(nb_atoms-1)/2))))
for step=1:nb_steps
    for i=1:nb_atoms
        for j=i+1:nb_atoms
            if cell_mod.distance(atoms[step],cell,i,j) < 1.7
                bond_matrix[step,i+j] += 1
            end
        end
    end
end

function autocorrelation(x,N,i,M)
    C=zeros(N)
    for k=i:i+M
        for n=1:N-1
            C[n] += x[k]*x[k-n]
        end
    end
    return C/C[1]
end

steps=5000
nb_steps=30000
time_corr_tot=zeros(steps)
for i=1:Int(trunc((nb_atoms*(nb_atoms-1)/2)))
    print("progress:",i/*(nb_atoms*(nb_atoms-1)/2)*100,"%\n")
    time_corr_loc=autocorrelation(bond_matrix[:,i],steps,steps,nb_steps-steps)
    for j=1:steps
        time_corr_tot[j] += time_corr_loc[j]
    end
end
for i=1:steps
    time_corr_tot[i] = time_corr_tot[i]/(nb_atoms*(nb_atoms-1)/2)
end

file_check=open(string("/home/moogmt/time_corr_",current_volume,"-",T,".dat"),"w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",time_corr_tot[Int(i)],"\n"))
end
close(file_check)
end
