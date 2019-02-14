GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75

min_lag=1
max_lag=5001
d_lag=5
unit=0.005


T=3000
V=8.82


folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"2-run/PIV_1/TRAJEC_test.xyz")

folder_out=string(folder_in,"Data/")

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_structure=size(traj)[1]

file_energy=open(string(folder_in,"2-run/EKS_base"))
lines=readlines(file_energy)
close(file_energy)

energy=zeros(nb_structure)
for i=1:nb_structure
    energy[i] = parse(Float64,split(lines[i])[3])
end

file_matrix=open(string(folder_in,"2-run/PIV_1/FRAME_TO_FRAME.MATRIX"))
lines=readlines(file_matrix)
close(file_matrix)

dist_max=parse(Float64,split(lines[1])[2])
distances=zeros(nb_structure,nb_structure)
for i=1:nb_structure
    print("Progress: ",i/nb_structure*100,"%\n")
    for j=i+1:nb_structure
        distances[i,j]=parse(Float64,split(lines[i+1])[j])*dist_max
        distances[j,i]=parse(Float64,split(lines[i+1])[j])*dist_max
    end
end

file_out=open(string(folder_in,"2-run/PIV_1/test.dat"),"w")
for i=1:nb_structure
    print("Progress: ",i/nb_structure*100,"%\n")
    for j=1:nb_structure
        write(file_out,string(distances[i,j]," ",abs(energy[i]-energy[j])/32,"\n"))
    end
end
close(file_out)
