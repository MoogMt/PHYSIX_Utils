GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Counts the number of bonds (and averages over time windows)

using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using cpmd
using markov


# Folder for data
folder_base="/media/moogmt/Stock/Gabriele/FRAME2FRAME/"

file=open(string(folder_base,"results_random.dat"))
lines=readlines(file)
close(file)

nb_structure=size(lines)[1]

index=zeros(nb_structure)

names_random=Array{AbstractString}(undef,nb_structure)
energy_ini=zeros(nb_structure)
energy_final=zeros(nb_structure)
for i=1:nb_structure
    names_random[i]        = split(lines[i])[1]
    energy_final[i] = parse(Float64,split(lines[i])[2])
    energy_ini[i]   = parse(Float64,split(lines[i])[6])
end

file=open(string(folder_base,"selection"))
lines=readlines(file)
close(file)
nb_structure_final=size(lines)[1]

names_final=Array{AbstractString}(undef,nb_structure_final)
for i=1:nb_structure_final
    line=split(lines[i])
    names_final[i] = line[3]
end

energy_ini_selection=zeros(nb_structure_final)
energy_final_selection=zeros(nb_structure_final)
for structure_selection=1:nb_structure_final
    for structure_random=1:nb_structure
        if names_random[structure_random] == names_final[structure_selection]
            energy_final_selection[structure_selection] = energy_final[structure_random]
            energy_ini_selection[structure_selection] = energy_ini[structure_random]
            break
        end
    end
end

# Determinig min
min_index=0
min_value=60
for i=1:nb_structure_final
    if min_value > energy_ini_selection[i]
        global min_value = energy_ini_selection[i]
        global min_index=i
    end
end
for i=1:nb_structure_final
    energy_ini_selection[i] -= min_value
end

file=open(string(folder_base,"FRAME_TO_FRAME.MATRIX_0"))
line=readline(file)
max_distance=parse(Float64,split(line)[2])
for i=1:min_index
    global line=readline(file)
end
line=split(line)
distances_PIV=zeros(nb_structure_final)
for i=1:nb_structure_final
    distances_PIV[i] = parse(Float64,line[i])*max_distance
end
close(file)

file_out=open(string(folder_base,"test_PIVvsE_ini.dat"),"w")
for i=1:nb_structure_final
    write(file_out,string(distances_PIV[i]," ",energy_ini_selection[i],"\n"))
end
close(file_out)


#===============================================================================#
#===============================================================================#

GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Gabriele/FRAME2FRAME_2/"

N_structure=[1,5,20]

for N in N_structure

    file=open(string(folder_base,"results_lowE_",N,".dat"))
    lines=readlines(file)
    close(file)

    nb_structure=size(lines)[1]
    nb_atoms=4

    index=zeros(nb_structure)

    names_random=Array{AbstractString}(undef,nb_structure)
    energy_ini=zeros(nb_structure)
    energy_final=zeros(nb_structure)
    for i=1:nb_structure
        names_random[i] = split(lines[i])[1]
        energy_final[i] = parse(Float64,split(lines[i])[2])/nb_atoms
        energy_ini[i]   = parse(Float64,split(lines[i])[6])/nb_atoms
    end

    # Determinig min
    min_index=0
    min_value=100000
    for i=1:nb_structure
        if min_value > energy_ini[i]
            min_value = energy_ini[i]
            min_index=i
        end
    end

    file=open(string(folder_base,"FRAME_TO_FRAME_",N,".MATRIX"))
    line=readline(file)
    max_distance=parse(Float64,split(line)[2])
    for i=1:min_index
        line=readline(file)
    end
    line=split(line)
    distances_PIV=zeros(nb_structure)
    for i=1:nb_structure
        distances_PIV[i] = parse(Float64,line[i])*max_distance
    end
    close(file)

    file_out=open(string(folder_base,"test_PIVvsE-",N,".dat"),"w")
    for i=1:nb_structure
        write(file_out,string(distances_PIV[i]," ",energy_ini[i],"\n"))
    end
    close(file_out)

    file_out=open(string(folder_base,"test_PIVvsDE-",N,".dat"),"w")
    for i=1:nb_structure
        write(file_out,string(distances_PIV[i]," ",energy_ini[i]-min_value,"\n"))
    end
    close(file_out)

end


#===============================================================================#
#===============================================================================#


GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"clustering.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Gabriele/FRAME2FRAME_2/"

N_structure=[1,5,20]

for N in N_structure

    file=open(string(folder_base,"results_lowE_",N,".dat"))
    lines=readlines(file)
    close(file)

    nb_structure=size(lines)[1]
    nb_atoms=4

    index=zeros(nb_structure)

    names_random=Array{AbstractString}(undef,nb_structure)
    phase_name=Array{AbstractString}(undef,nb_structure)
    energy_ini=zeros(nb_structure)
    energy_final=zeros(nb_structure)
    for i=1:nb_structure
        line=split(lines[i])
        names_random[i] = line[1]
        phase_name[i]   = line[4]
        energy_final[i] = parse(Float64,line[2])/nb_atoms
        energy_ini[i]   = parse(Float64,line[6])/nb_atoms
    end

    # Determinig min
    min_index=0
    min_value=100000
    for i=1:nb_structure
        if min_value > energy_ini[i]
            min_value = energy_ini[i]
            min_index=i
        end
    end

    file=open(string(folder_base,"FRAME_TO_FRAME_",N,".MATRIX"))
    line=readline(file)
    max_distance=parse(Float64,split(line)[2])
    for i=1:min_index
        line=readline(file)
    end
    line=split(line)
    distances_PIV=zeros(nb_structure)
    for i=1:nb_structure
        distances_PIV[i] = parse(Float64,line[i])*max_distance
    end
    close(file)

    file_out=open(string(folder_base,"test_PIVvsE-",N,".dat"),"w")
    for i=1:nb_structure
        write(file_out,string(distances_PIV[i]," ",energy_ini[i],"\n"))
    end
    close(file_out)

    file_out=open(string(folder_base,"test_PIVvsDE-",N,".dat"),"w")
    for i=1:nb_structure
        write(file_out,string(distances_PIV[i]," ",energy_ini[i]-min_value,"\n"))
    end
    close(file_out)


    indexs_ranking=clustering.simpleSequence(nb_structure)
    for i=1:nb_structure
        for j=i+1:nb_structure
            if energy_ini[indexs_ranking[i]] > energy_ini[indexs_ranking[j]]
                stock=indexs_ranking[i]
                indexs_ranking[i]=indexs_ranking[j]
                indexs_ranking[j]=stock
            end
        end
    end

    rho=zeros(nb_structure)
    for i=1:nb_structure
        rho[i] = energy_ini[indexs_ranking[i]]-energy_ini[indexs_ranking[1]]
    end

    file=open(string(folder_base,"FRAME_TO_FRAME_",N,".MATRIX"))
    line=readline(file)
    max_distance=parse(Float64,split(line)[2])
    close(file)

    # TO DO : Implement Burning
    dc_burn=0.5
    burn=zeros(nb_structure)
    for target_structure=1:nb_structure
        file=open(string(folder_base,"FRAME_TO_FRAME_",N,".MATRIX"))
        line=readline(file)
        for i=1:indexs_ranking[target_structure]
            line=readline(file)
        end
        line=split(line)
        for second_structure=1:nb_structure
            if second_structure == target_structure
                continue
            end
            if  dc_burn > parse(Float64,line[indexs_ranking[second_structure]])*max_distance && energy_ini[indexs_ranking[second_structure]] <  energy_ini[indexs_ranking[target_structure]]
                burn[indexs_ranking[target_structure]] = 1
            end
        end
    end

    delta=ones(nb_structure)*max_distance
    nn_neigh=zeros(Int,nb_structure)
    for target_structure=1:nb_structure
        file=open(string(folder_base,"FRAME_TO_FRAME_",N,".MATRIX"))
        line=readline(file)
        for i=1:indexs_ranking[target_structure]
            line=readline(file)
        end
        line=split(line)
        for second_structure=1:nb_structure
            if second_structure == target_structure
                continue
            end
            dist=parse(Float64,line[indexs_ranking[second_structure]])*max_distance
            if dist < delta[indexs_ranking[target_structure]] && energy_ini[indexs_ranking[second_structure]] < energy_ini[indexs_ranking[target_structure]] && burn[second_structure] < 1
                delta[indexs_ranking[target_structure]] = dist
                nn_neigh[indexs_ranking[target_structure]] = indexs_ranking[second_structure]
            end
        end
        close(file)
    end

    delta[1]=max_distance

    file_out=open(string(folder_base,"decision_diagram-",N,".dat"),"w")
    for i=1:nb_structure
        if burn[i] < 1
            write(file_out,string(rho[i]," ",delta[i],"\n"))
        end
    end
    close(file_out)

    min_delta=1.2
    max_rho=4

    n_cluster=0
    cl=ones(Int,nb_structure)*(-1)
    icl=[]
    # Determine the cluster centers
    for i=1:nb_structure
        if rho[i] < max_rho && delta[i] > min_delta && burn[i] < 1
    		n_cluster += 1
            cl[indexs_ranking[i]] = n_cluster
            icl=push!(icl,indexs_ranking[i])
        end
    end

    for i=1:nb_structure
        if cl[indexs_ranking[i]] == -1
            cl[indexs_ranking[i]] = cl[ indexs_ranking[nn_neigh[i]]  ]
        end
    end
    for i=1:n_cluster
        file=open(string(folder_base,"cluster",i,"-",N,".dat"),"w")
        for j=1:nb_structure
            if i == cl[j]
                write(file,string(distances_PIV[indexs_ranking[j]]," ",rho[indexs_ranking[j]]," ",delta[indexs_ranking[j]],"\n"))
            end
        end
        close(file)
    end


    file_out2=open(string(folder_base,"PhaseNames-",N,".dat"),"w")
    for cluster=1:n_cluster
        file_out=open(string(folder_base,"check_cluster-",N,"-",cluster,".dat"),"w")
        name_check=phase_name[icl[cluster]]
        check_rho=10
        check_delta=0
        for i=1:nb_structure
            if name_check == phase_name[indexs_ranking[i]]
                if rho[indexs_ranking[i]] < check_rho
                    check_rho = rho[indexs_ranking[i]]
                    check_delta = delta[indexs_ranking[i]]
                end
            end
        end
        write(file_out,string(check_rho," ",check_delta,"\n"))
        write(file_out2,string(cluster," ",name_check,"\n"))
        close(file_out)
    end
    close(file_out2)

end
