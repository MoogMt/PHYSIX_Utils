GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

cut_off_states=0.1

min_lag=1
max_lag=5001
d_lag=5
unit=0.005

# for V in Volumes
#     for T in Temperatures


folder_in=string(folder_base,V,"/",T,"K/")
folder_out=string(folder_in,"Data/")

# if  ! isfile(string(folder_out,"TransitionsMatrix-C.dat"))
#     continue
# end

transition_matrix=readTransitionMatrix(string(folder_out,"TransitionsMatrix-C.dat"))
file=open(string(folder_out,"test_life.dat"),"w")
for step=1:size(transition_matrix)[3]
        write(file,string(step," ",1/(1-transition_matrix[1,1,step])*unit*d_lag,"\n"))
end
close(file)

#     end
# end
