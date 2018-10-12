include("contactmatrix.jl")

folder=string("/media/moogmt/Stock/MoS2/Relax/PBE/Mo2S4/")
min_energy=0
max_cluster=1000

Ry2eV=13.6056980659

file_candidates=open(string(folder,"candidates.xyz"),"w")
file_energy_candidates=open(string(folder,"candidates_energy.dat"),"w")

for i=1:max_cluster
    folder_local=string(folder,"/",i,"_center/")
    file=string("cluster_candidate.xyz")
    if isfile( string(folder_local,file) ) && (! isfile(string(folder_local,"CRASH")))
        # ( QE prints another structure as final step for some reason )

        print("Dealing with ",i," cluster\n")
        # Getting the last step of the relax
        traj = filexyz.readFastFile(string(folder_local,file))
        if size(traj)[1] == 0
            continue
        end
        cluster_candidate=traj[size(traj)[1]-1]

        # Energy
        # Reading file
        file_energy=open(string(folder_local,"energy"))
        lines=readlines(file_energy);
        close(file_energy)
        # Getting next to last step
        energy=parse(Float64,lines[size(lines)[1]-1])

        # Keeping in mind what is the min energy
        if min_energy > energy
            min_energy = energy
        end
        # Writting structure to file
        write(file_candidates,string(size(cluster_candidate.names)[1],"\n"))
        write(file_candidates,string(" cluster: ",i,"\n"))
        for atom=1:size(cluster_candidate.names)[1]
            write(file_candidates,string(cluster_candidate.names[atom]," "))
            for j=1:3
                write(file_candidates,string(cluster_candidate.positions[atom,j]," "))
            end
            write(file_candidates,string("\n"))
        end
        # Writting energy to file
        write(file_energy_candidates,string(i," ",energy*Ry2eV,"\n"))
    end
end
close(file_candidates)
close(file_energy_candidates)

file_min_energy=open(string(folder,"min_energy.dat"),"w")
write(file_min_energy,string(min_energy*Ry2eV,"\n"))
close(file_min_energy)

min_energy=min_energy*Ry2eV

# Reading the files
clusters = filexyz.readFastFile(string(folder,"candidates.xyz"))
file_cluster=open(string(folder,"candidates_energy.dat"))
lines=readlines(file_cluster);
close(file_cluster)

# Reconstructing structures
nb_structure=size(lines)[1]
energy=zeros(nb_structure)
for i=1:nb_structure
    energy[i]= parse(Float64,split(lines[i])[2])
end


# Sorting by increasing energy
for i=1:nb_structure
    for j=1:nb_structure
        if energy[i] < energy[j]
            stock=energy[i]
            energy[i]=energy[j]
            energy[j]=stock
            stock_cluster=clusters[i]
            clusters[i]=clusters[j]
            clusters[j]=stock_cluster
        end
    end
end

# Check
debug=false
if debug
    sorted=open(string(folder,"sorted_energy.dat"),"w")
    for i=1:nb_structure
        write(sorted,string(energy[i],"\n"))
    end
    close(sorted)
end


#=================================================================================#

folder=string("/media/moogmt/Stock/MoS2/Relax/PBE/Mo3S6/")
min_energy=0
max_cluster=1000

file_candidates_3=open(string(folder,"candidates.xyz"),"w")
file_energy_candidates_3=open(string(folder,"candidates_energy.dat"),"w")

for cl=1:4
    for i=1:max_cluster
        folder_local=string(folder,"CL",cl,"/",i,"_center/")
        file=string("cluster_candidate.xyz")
        if isfile( string(folder_local,file) )
            # ( QE prints another structure as final step for some reason )

            # Getting the last step of the relax
            traj = filexyz.readFastFile(string(folder_local,file))
            cluster_candidate=traj[size(traj)[1]-1]

            # Energy
            # Reading file
            file_energy=open(string(folder_local,"energy"))
            lines=readlines(file_energy);
            close(file_energy)
            # Getting next to last step
            energy=parse(Float64,lines[size(lines)[1]-1])

            # Keeping in mind what is the min energy
            if min_energy > energy
                min_energy = energy
            end
            # Writting structure to file
            write(file_candidates_3,string(size(cluster_candidate.names)[1],"\n"))
            write(file_candidates_3,string("CL: ",cl," cluster: ",i,"\n"))
            for atom=1:size(cluster_candidate.names)[1]
                write(file_candidates_3,string(cluster_candidate.names[atom]," "))
                for j=1:3
                    write(file_candidates_3,string(cluster_candidate.positions[atom,j]," "))
                end
                write(file_candidates_3,string("\n"))
            end
            # Writting energy to file
            write(file_energy_candidates_3,string(cl," ",i," ",energy*Ry2eV,"\n"))
        end
    end
end
close(file_candidates_3)
close(file_energy_candidates_3)

file_min_energy=open(string(folder,"min_energy.dat"),"w")
write(file_min_energy,string(min_energy*Ry2eV,"\n"))
close(file_min_energy)

#==============================================================================#

folder=string("/media/moogmt/Stock/MoS2/Relax/PBE/Mo5S10/")
min_energy=0
max_cluster=1000

Ry2ev=13.6056980659

file_candidates=open(string(folder,"candidates.xyz"),"w")
file_energy_candidates=open(string(folder,"candidates_energy.dat"),"w")

for i=1:max_cluster
    folder_local=string(folder,"/",i,"_center/")
    file=string("cluster_candidate.xyz")
    if isfile( string(folder_local,file) ) && (! isfile(string(folder_local,"CRASH")))
        # ( QE prints another structure as final step for some reason )

        print("Dealing with ",i," cluster\n")
        # Getting the last step of the relax
        traj = filexyz.readFastFile(string(folder_local,file))
        if size(traj)[1] == 0
            continue
        end
        cluster_candidate=traj[size(traj)[1]-1]

        # Energy
        # Reading file
        file_energy=open(string(folder_local,"energy"))
        lines=readlines(file_energy);
        close(file_energy)
        # Getting next to last step
        energy=parse(Float64,lines[size(lines)[1]-1])

        # Keeping in mind what is the min energy
        if min_energy > energy
            min_energy = energy
        end
        # Writting structure to file
        write(file_candidates,string(size(cluster_candidate.names)[1],"\n"))
        write(file_candidates,string(" cluster: ",i,"\n"))
        for atom=1:size(cluster_candidate.names)[1]
            write(file_candidates,string(cluster_candidate.names[atom]," "))
            for j=1:3
                write(file_candidates,string(cluster_candidate.positions[atom,j]," "))
            end
            write(file_candidates,string("\n"))
        end
        # Writting energy to file
        write(file_energy_candidates,string(i," ",energy*Ry2eV,"\n"))
    end
end
close(file_candidates)
close(file_energy_candidates)

file_min_energy=open(string(folder,"min_energy.dat"),"w")
write(file_min_energy,string(min_energy*Ry2eV,"\n"))
close(file_min_energy)

#===============================================================================#
