GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

V=9.8
T=3000

nbC=32
nbO=nbC*2

cut_off_bond = 1.75
max_neigh=5

min_lag=1
max_lag=5001
d_lag=5
unit=0.005

# for V in Volumes
#     for T in Temperatures

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]

states, state_matrix, percent = assignDataToStates( data, nb_types, type_atoms )

case_states=size(states[1])[2]
nb_states=size(states[1])[1]
normC=4
normO=2
oks_states=zeros(nb_states)
for state=1:nb_states
    part1=4
    part2=0
    for i=1:case_states
        if states[1][state,i] > 0
            part1 -= 1
            if i > max_neigh # O
                part2 += (normO-states[1][state,i])
            else # C
                part2 += (normC-states[1][state,i])
            end
        end
    end
    if part1 == part2
        oks_states[state] = 1
    end
end

nb_states_exc=0
count_states_all=0
for i=1:nb_states
    if oks_states[i] != 1
        global nb_states_exc += percent[1][i]
    end
    global count_states_all+= percent[1][i]
end
percent_excited=nb_states_exc/count_states_all

coordinancesC=zeros(max_neigh)
for state=1:nb_states
    nb=0
    for i=1:case_states
        if states[1][state,i] > 0
            nb += 1
        end
    end
    coordinancesC[nb] += percent[1][state]
end
for i=1:max_neigh
    coordinancesC[i] /= count_states_all
end
coordinancesC *= 100

percent_CO2=0
stateC=[0,0,0,0,0,1,1,0,0,0]
for state=1:nb_states
    found=true
    for i=1:case_states
        if stateC[i] != states[1][state,i]
            found=false
        end
    end
    if found
        global percent_CO2=percent[1][state]/count_states_all*100
        break
    end
end


stateC=[0,0,0,0,0,2,2,1,0,0]
nb_state=0
for state=1:nb_states
    found=true
    for i=1:case_states
        if stateC[i] != states[1][state,i]
            found=false
        end
    end
    if found
        global nb_state=state
        break
    end
end


for carbon=1:nbC
    for step=1:nb_steps
        for
        end
    end
end



stateO=[3,3,0,0,0,0,0,0,0,0]
stateC=[0,0,0,0,0,2,2,0,0,0]

O_nb=0
for i=1:size(states)[1]
    dO=0
    for j=1:size(states)[2]
        dO += (states[i,j]-stateO[j])^2
    end
    if dO == 0
        O_nb = i
    end
end

taumax=100
dtau=1
file_out=open(string(folder_out,"Test.dat"),"w")
for step=1:size(state_matrix)[2]
    print("Progress: ",step/size(state_matrix)[2]*100,"%\n")
    count=0
    for oxygen=1:nbC
        if state_matrix[nbC+oxygen,step] == O_nb
            for oxygen2=oxygen+1:nbO
                if state_matrix[nbC+oxygen2,step] == O_nb
                    for carbon=1:nbC
                        if cell_mod.distance(traj[step],cell,carbon,nbC+oxygen) < cut_off_bond && cell_mod.distance(traj[step],cell,carbon,nbC+oxygen2) < cut_off_bond
                            for carbon2=carbon+1:nbC
                                if cell_mod.distance(traj[step],cell,carbon2,nbC+oxygen) < cut_off_bond && cell_mod.distance(traj[step],cell,carbon2,nbC+oxygen2) < cut_off_bond
                                    count += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    write(file_out,string(step," ",count,"\n"))
end
close(file_out)
