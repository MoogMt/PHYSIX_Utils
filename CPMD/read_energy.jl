GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod

avg_times=zeros(3)
avg_var=zeros(3)

#==============================================================================#
nb_run_2=4
for i=1:nb_run_2
    #
    folder_base=string("/media/moogmt/Stock/Mathieu/MoS2/MetaExplore/CPMD/Mo2S4/run_",i,"/")
    file_in=open(string(folder_base,"ENERGIES"))
    lines=readlines(file_in)
    close(file_in)
    #
    nb_lines=size(lines)[1]
    avg_time=0
    var_time=0
    for i=1:nb_lines
        max_col=size(split(lines[i]))[1]
        time_ = parse(Float64,split(lines[i])[max_col])
        avg_time += time_
        var_time += time_*time_
    end
    avg_time=avg_time/nb_lines
    var_time=sqrt( var_time/nb_lines - avg_time*avg_time )
    avg_times[1] += avg_time
    avg_var[1] += var_time
end
avg_times[1]=avg_times[1]/nb_run_2
avg_var[1] = avg_var[1]/nb_run_2
#==============================================================================#
nb_run_3=4
for i=1:nb_run_3
    folder_base=string("/media/moogmt/Stock/Mathieu/MoS2/MetaExplore/CPMD/Mo3S6/run_",1,"/")
    file_in=open(string(folder_base,"ENERGIES"))
    lines=readlines(file_in)
    close(file_in)
    nb_lines=size(lines)[1]
    avg_time=0
    var_time=0
    for i=1:nb_lines
        max_col=size(split(lines[i]))[1]
        time_ = parse(Float64,split(lines[i])[max_col])
        avg_time += time_
        var_time += time_*time_
    end
    avg_time=avg_time/nb_lines
    var_time=sqrt( var_time/nb_lines - avg_time*avg_time )
    avg_times[2] += avg_time
    avg_var[2] += var_time
end
avg_times[2]=avg_times[2]/nb_run_3
avg_var[2] = avg_var[2]/nb_run_3
#==============================================================================#
nb_run_4=10
for i=1:nb_run_4
    folder_base=string("/media/moogmt/Stock/Mathieu/MoS2/MetaExplore/CPMD/Mo4S8/run_",i,"/")
    file_in=open(string(folder_base,"ENERGIES"))
    lines=readlines(file_in)
    close(file_in)
    nb_lines=size(lines)[1]
    avg_time=0
    var_time=0
    for i=1:nb_lines
        max_col=size(split(lines[i]))[1]
        time_ = parse(Float64,split(lines[i])[max_col])
        avg_time += time_
        var_time += time_*time_
    end
    avg_time=avg_time/nb_lines
    var_time=sqrt( var_time/nb_lines - avg_time*avg_time )
    avg_times[3] += avg_time
    avg_var[3] += var_time
end
avg_times[3]=avg_times[3]/nb_run_3
avg_var[3] = avg_var[3]/nb_run_3
#==============================================================================#

file_out=open(string("/media/moogmt/Stock/Mathieu/MoS2/MetaExplore/CPMD/time.dat"),"w")
for i=1:3
    write(file_out,string(i+1," ",avg_times[i]," ",avg_var[i],"\n"))
end
close(file_out)
