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
    N=Int(trunc(M*frac)-1)
    C=zeros(N)
    for n=1:N
        count=0
        for m=1:M-n+1
            C[n] += x[m]*x[m+(n-1)]
            count +=1
        end
        C[n] = C[n]/count
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

for carbon in [1]
    for oxygen in [35-32]
        bond_vector=zeros(nb_steps)
        mean_value=0
        for step=1:nb_steps
            print("Bond function - C:",carbon," O:",oxygen,"\n")
            if cell_mod.distance(atoms[step],cell,carbon,nbC+oxygen) < cut_off
                bond_vector[step] = 1
                mean_value += 1
            end
        end
        mean_value /= nb_steps
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


#
# only_bond=false
#
# if only_bond
#
#     atoms1=[]
#     atoms2=[]
#     for i=1:32
#         for j=33:nb_atoms
#             if cell_mod.distance(atoms[1],cell,i,j) < cutoff
#                 push!(atoms1,i)
#                 push!(atoms2,j)
#             end
#         end
#     end
#
#
#     for i=1:size(atoms1)[1]
#         bond_matrix=zeros(nb_steps)
#         for step=1:nb_steps
#             if cell_mod.distance(atoms[step],cell,atoms1[i],atoms2[i]) < cutoff
#                 bond_matrix[step] = 1
#             end
#         end
#         corr=autocorrelation2(bond_matrix,frac)
#         for i=1:steps-1
#             time_corr[i] += corr[i]
#         end
#     end
#
#     movingOs=[]
#     targetCs=[]
#     file_dist=open(string("/home/moogmt/dist.dat"),"w")
#     for step=1:nb_steps-1
#         write(file_dist,string(step*stride*unit," "))
#         for i=1:size(atoms1)[1]
#             write(file_dist,string(cell_mod.distance(atoms[step],cell,atoms1[i],atoms2[i])-cell_mod.distance(atoms[1],cell,atoms1[i],atoms2[i])," "))
#             if cell_mod.distance(atoms[step],cell,atoms1[i],atoms2[i])-cell_mod.distance(atoms[1],cell,atoms1[i],atoms2[i]) > 1
#                 print("check ",atoms1[i]," ",atoms2[i],"\n")
#             end
#         end
#         write(file_dist,"\n")
#     end
#     close(file_dist)
#
#     file_atoms=open(string("/home/moogmt/check_atoms.dat"),"w")
#     print("CHECK:", size(movingOs)[1],"\n")
#     for w=1:size(movingOs)[1]
#         write(file_atoms,string(movingOs[w]," ",targetCs[w],"\n"))
#     end
#     close(file_atoms)
#
# else
#
#     for atom1=1:32
#         print("progress:",atom1/32*100,"%\n")
#         for atom2=33:nb_atoms
#             bond_matrix=zeros(nb_steps)
#             for step=1:nb_steps
#                 if cell_mod.distance(atoms[step],cell,atom1,atom2) < cutoff
#                     bond_matrix[step] = 1
#                 end
#             end
#             corr=autocorrelation2(bond_matrix,frac)
#             for i=1:steps-1
#                 time_corr[i] += corr[i]
#             end
#         end
#     end
#
# end

# Normalization
# value=time_corr[1]
# time_corr /= value
#
# value=corr_stock[1]
# corr_stock /= value
#
# print("atom1:",index1," atom2: ",index2,"\n")
#
# file_check=open(string("/home/moogmt/time_long.dat"),"w")
# for i=1:steps-1
#     write(file_check,string(i*unit*stride," ",corr_stock[Int(i)],"\n"))
# end
# close(file_check)

# end
#end
#
# file_dist=open(string("/home/moogmt/dist.dat"),"w")
# for step=1:nb_steps-1
#     write(file_dist,string(step*stride*unit," ",cell_mod.distance(atoms[step],cell,8,42),"\n"))
# end
# close(file_dist)
