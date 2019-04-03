GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

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
    if min_value > energy_final_selection[i]
        global min_value = energy_final_selection[i]
        global min_index=i
    end
end
for i=1:nb_structure_final
    energy_final_selection[i] -= min_value
end

file=open(string(folder_base,"FRAME_TO_FRAME.MATRIX_100"))
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

file_out=open(string(folder_base,"test_PIVvsE_final.dat"),"w")
for i=1:nb_structure_final
    write(file_out,string(distances_PIV[i]," ",energy_final_selection[i],"\n"))
end
close(file_out)

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
    if min_value > energy_final_selection[i]
        global min_value = energy_final_selection[i]
        global min_index=i
    end
end
for i=1:nb_structure_final
    energy_final_selection[i] -= min_value
end
# Determinig min
min_value=60
for i=1:nb_structure_final
    if min_value > energy_final_selection[i]
        global min_value = energy_ini_selection[i]
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

file_out=open(string(folder_base,"test_PIVvsE_ini_RefSkewed.dat"),"w")
for i=1:nb_structure_final
    write(file_out,string(distances_PIV[i]," ",energy_ini_selection[i],"\n"))
end
close(file_out)
