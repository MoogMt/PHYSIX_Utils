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

for T in Temperatures
    file_coordinanceC=open(string(folder_base,"coordinancesC.dat"),"w")
    file_coordinanceO=open(string(folder_base,"coordinancesC.dat"),"w")
    file_coordinanceC=open(string(folder_base,"coordinancesC.dat"),"w")
    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"TRAJEC_wrapped.xyz")
        folder_out=string(folder_in,"Data/")

        if ! isfile(file)
            continue
        end

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
        nb_state_dimer=0
        for state=1:nb_states
            found=true
            for i=1:case_states
                if stateC[i] != states[1][state,i]
                    found=false
                end
            end
            if found
                global nb_state_dimer=state
                break
            end
        end


        dt=25
        time_stock=[]
        count_exchange=0
        carbon_stock=[]
        for carbon=1:nbC
            record=zeros(nb_steps)
            for step=1:nb_steps
                found_one=false
                if state_matrix[1][carbon,step] == nb_state_dimer
                    # Check the exchange in two step:
                    max_neigh2=zeros(2)
                    count_bond_main=0
                    for atom=1:nb_atoms
                        if atom == carbon
                            continue
                        end
                        step_before=step-dt
                        if step_before <= 0
                            step_before=1
                        end
                        step_after = step+dt
                        if step_after > nb_steps
                            step_after=nb_steps
                        end
                        if cell_mod.distance(traj[step_before],cell,carbon,atom) < cut_off_bond && cell_mod.distance(traj[step_after],cell,carbon,atom) < cut_off_bond
                            count_bond_main += 1
                        end
                    end
                    if count_bond_main < 2 && record[step] < 1
                        for i=1:dt
                            if step+dt < nb_steps
                                record[step+i] += 1
                            else
                                break
                            end
                        end
                        global count_exchange += 1
                        push!(time_stock,step)
                        push!(carbon_stock,carbon)
                        found_one=true
                    end
                end
                if found_one
                    if step+dt < nb_steps
                        step += dt
                    else
                        break
                    end
                end
            end
        end

        count_exchange2=0
        carbon_exchange=zeros(nbC)
        time_stock2=[]
        carbon_stock2=[]
        check=zeros(size(carbon_stock)[1])
        for i=1:size(carbon_stock)[1]-1
            print("size:",size(check)[1],"\n")
            found=false
            marker=0
            for j=i+1:size(carbon_stock)[1]
                if check[j] == true
                    continue
                end
                if time_stock[i] == time_stock[j]
                    check[i] = 1
                    check[j] = 1
                    push!(time_stock2,time_stock[i])
                    push!(carbon_stock2,carbon_stock[i])
                    push!(time_stock2,time_stock[j])
                    push!(carbon_stock2,carbon_stock[j])
                    found=true
                    marker=j
                    break
                end
            end
            if found
                global count_exchange += 1
                if carbon_exchange[j] < 1
                    carbon_exchange[j] += 1
                end
                if carbon_exchange[i] < 1
                    carbon_exchange[i] += 1
                end
            end
        end
        nb_carbons=sum(carbon_exchange)
        frac_carbon=nb_carbons/nbC*100

    end
end
