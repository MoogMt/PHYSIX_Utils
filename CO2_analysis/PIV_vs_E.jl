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
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")

print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

n_structure=4001

traj=traj[1:n_structure]

file_energy=open(string(folder_in,"1-run/PIV_test/eks"))
lines=readlines(file_energy)
close(file_energy)

energy=zeros(n_structure)
for i=1:n_structure
    energy[i] = parse(Float64,split(lines[i])[3])
end

file_matrix=open(string(folder_in,"1-run/PIV_test/FRAME_TO_FRAME.MATRIX"))
lines=readlines(file_matrix)
close(file_matrix)


n_test=1000
distances=zeros(n_test,n_structure)
for i=1:n_test
    print("Progress: ",i/n_test*100,"%\n")
    for j=1:n_structure
        distances[i,j]=parse(Float64,split(lines[i+1])[j])
    end
end

file_out=open(string(folder_in,"1-run/PIV_test/test.dat"),"w")
for i=1:n_test
    print("Progress: ",i/n_test*100,"%\n")
    for j=1:n_structure
        if abs(energy[i]-energy[j])*13.6 < 2 && distances[i,j] < 0.3
            write(file_out,string(distances[i,j]," ",abs(energy[i]-energy[j]),"\n"))
        end
    end
end
close(file_out)
