include("contactmatrix.jl")

# Thermodynamical values
Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]
cut_offs=[1.6,1.75,1.8]

# Current Volume and Temperature
current_volume=parse(Float64,Volumes[1])
current_temperature=Temperature[4]
folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",current_volume,"/3000K/")
file=string(folder,"TRAJEC_wrapped.xyz")

# Time values
unit=0.0005*5
stride=1

atoms = filexyz.read(file,stride)
cell=cell_mod.Cell_param(current_volume,current_volume,current_volume)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

cutoff=1.8
atom=1
bond_matrix=zeros(nb_steps,nb_atoms)
for i=1:nb_steps
    for j=1:nb_atoms
        if j != atom
            if cell_mod.distance(atoms[i],cell,atom,j) < co
                bond_matrix[i,j] += 1
                print(bond_matrix[i,j])
            end
        end
    end
end


steps=Int(trunc(nb_steps*0.01))
time_corr_tot=zeros(steps)
for i=32:nb_atoms
    print("atom:",i,"\n")
    time_corr_loc=zeros(steps)
    for t0=1:steps
        for t=1:nb_steps-1000-t0
            print("check:",bond_matrix[t,i]*bond_matrix[t+t0,i],"\n")
            time_corr_loc[t0] += bond_matrix[t,i]*bond_matrix[t+t0,i]
        end
        time_corr_loc[t0] /= (steps-t0)
    end
    time_corr_tot += time_corr_loc
end
time_corr_tot /= nb_atoms

file_check=open("/home/moogmt/time_corr.dat","w")
for i=1:steps
    write(file_check,string(i*unit*stride," ",time_corr_tot[Int(i)],"\n"))
end
close(file_check)
