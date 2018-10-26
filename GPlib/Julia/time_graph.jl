include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]
Cut_Off=[1.75]

# Time values
unit=0.005
frac=0.4

# Number of atoms
nbC=32
nbO=nbC*2

# for cut_off in Cut_Off
#     for T in Temperatures
#         for V in Volumes

V=9.8
T=3000
cut_off=1.75


folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

if ! isfile(folder_in,"TRAJEC_wrapped.xyz")
    # continue
    print("Check\n")
end

folder_out=string(folder_in,"Data/")

atoms = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_atoms=size(atoms[1].names)[1]
nb_steps=size(atoms)[1]
steps=Int(trunc(nb_steps*frac))

time_corr_gen=zeros(steps)
time_corr_bonds=zeros(nbC,nbO,steps)

for carbon=1:nbC
    for oxygen=1:nbO
        bond_vector=zeros(nb_steps)
        mean_value=0
        for step=1:nb_steps
            print("Bond function - C:",carbon," O:",oxygen," Progress:",step/nb_steps*100,"%\n")
            if cell_mod.distance(atoms[step],cell,carbon,nbC+oxygen) < cut_off
                bond_vector[step] = 1
                mean_value += 1
            end
        end
        mean_value /= nb_steps
        # If there is no bond except for flickering, we just don't care
        if sum(bond_vector) < 10
            continue
        end
        # Autocorrelation
        for i=0:steps-1
            print("Autocorrelation - C:",carbon," O:",oxygen," Progress: ",i/(steps-1)*100,"%\n")
            for j=1:nb_steps-i-1
                time_corr_bonds[carbon,oxygen,i+1]+=bond_vector[j]*bond_vector[j+i]
            end
            time_corr_bonds[carbon,oxygen,i+1] /= nb_steps-i
            time_corr_gen += time_corr_bonds[carbon,oxygen,i+1]
        end
    end
end
time_corr_gen /= nbC*nbO


file_check=open(string(folder_out,"bond_correlation_general.dat"),"w")
for i=1:steps-1
    write(file_check,string(i*unit," ",time_corr_gen[i]/nb_steps,"\n"))
end
close(file_check)

file_check=open(string(folder_out,"bond_correlation_individual.dat"),"w")
for i=1:steps-1
    print("Writting to file - Progress: ",i/(steps-1)*100,"%\n")
    write(file_check,string(i*unit," "))
    for carbon=1:nbC
        for oxygen=1:nbO
            write(file_check,string(time_corr_bonds[carbon,oxygen,i]," "))
        end
    end
    write(file_check,string("\n"))
end
close(file_check)
