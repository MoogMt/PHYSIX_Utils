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

# Counting the number of structure below a threshold of energy (eV)
energy_threshold=2
count=0
for i=1:nb_structure
    if abs(energy[i]-min_energy) < 2
        count += 1
    end
end
# From that point we care only about the structures below the threshold
nb_structure = count

# Building PIV-like vector for each cluster
nb_atoms=size(clusters[1].names)[1]
n_mo=Int(nb_atoms/3)
size_mo=Int(n_mo*(n_mo-1)/2)
n_s=Int(nb_atoms/3*2)
size_s=Int(n_s*(n_s-1)/2)
size_mos=n_s*n_mo
size_total=Int(nb_atoms*(nb_atoms-1)/2)
matrix=zeros(nb_structure, size_total )
for cluster=1:nb_structure
    #================================================#
    # compute Mo-Mo part
    for i=1:n_mo-1
        for j=i+1:n_mo
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])
            end
            matrix[cluster,(i-1)*n_mo+j-1]=sqrt(dist)
        end
    end
    for i=1:size_mo-1
        for j=i+1:size_mo
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    # Compute Mo-S part
    for i=1:n_mo
        for j=1:n_s
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+1+(i-1)*n_s+(j-1)]=sqrt(dist)
        end
    end
    for i=size_mo+1:size_mo+size_mos
        for j=i+1:size_mo+size_mos+1
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    #Compute S-S part
    test=1
    for i=1:n_s-1
        for j=i+1:n_s
            dist=0
            for k=1:3
                dist=(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+size_mos+test]=sqrt(dist)
            test+=1
        end
    end
    for i=size_mo+size_mos+1:size_total-1
        for j=i+1:size_total
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
end

# Striking structure too close to those already existant
cut_off_distance=1.5 # Cut_off for topological distance, in Angstrom
strike=zeros(nb_structure)
for i=1:nb_structure-1
    for j=i+1:nb_structure
        dist=0
        for k=1:size_total
            dist += sqrt((matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k]))
        end
        if dist/size_total < cut_off_distance
            strike[j]=1
        end
    end
end

# Writting remaining structures
energy_final=open(string(folder,"energy_final.dat"),"w")
xyz_final=open(string(folder,"all_relax_cluster.xyz"),"w")
count=1
for i=1:nb_structure
    if strike[i] == 0
        write(energy_final,string(energy[i],"\n"))
        filexyz.write(xyz_final,clusters[i])
        file_xyz_loc=open(string(folder,count,"_cluster_relax.xyz"),"w")
        filexyz.write(file_xyz_loc,clusters[i])
        close(file_xyz_loc)
        count += 1
    end
end
close(energy_final)
close(xyz_final)

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
    energy[i]= parse(Float64,split(lines[i])[3])
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

# Counting the number of structure below a threshold of energy (eV)
energy_threshold=2
count=0
for i=1:nb_structure
    if abs(energy[i]-min_energy) < 2
        count += 1
    end
end
# From that point we care only about the structures below the threshold
nb_structure = count
print("COUNT: ",count,"\n")

# Building PIV-like vector for each cluster
nb_atoms=size(clusters[1].names)[1]
n_mo=Int(nb_atoms/3)
size_mo=Int(n_mo*(n_mo-1)/2)
n_s=Int(nb_atoms/3*2)
size_s=Int(n_s*(n_s-1)/2)
size_mos=n_s*n_mo
size_total=Int(nb_atoms*(nb_atoms-1)/2)
matrix=zeros(nb_structure, size_total )
for cluster=1:nb_structure
    #================================================#
    # compute Mo-Mo part
    for i=1:n_mo-1
        for j=i+1:n_mo
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])
            end
            matrix[cluster,(i-1)*n_mo+j-1]=sqrt(dist)
        end
    end
    for i=1:size_mo-1
        for j=i+1:size_mo
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    # Compute Mo-S part
    for i=1:n_mo
        for j=1:n_s
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+1+(i-1)*n_s+(j-1)]=sqrt(dist)
        end
    end
    for i=size_mo+1:size_mo+size_mos
        for j=i+1:size_mo+size_mos+1
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    #Compute S-S part
    test=1
    for i=1:n_s-1
        for j=i+1:n_s
            dist=0
            for k=1:3
                dist=(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+size_mos+test]=sqrt(dist)
            test+=1
        end
    end
    for i=size_mo+size_mos+1:size_total-1
        for j=i+1:size_total
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
end

# Striking structure too close to those already existant
cut_off_distance=1.5 # Cut_off for topological distance, in Angstrom
strike=zeros(nb_structure)
for i=1:nb_structure-1
    for j=i+1:nb_structure
        dist=0
        for k=1:size_total
            dist += sqrt((matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k]))
        end
        if dist/size_total < cut_off_distance
            strike[j]=1
        end
    end
end

# Writting remaining structures
energy_final=open(string(folder,"energy_final.dat"),"w")
xyz_final=open(string(folder,"all_relax_cluster.xyz"),"w")
count=1
for i=1:nb_structure
    if strike[i] == 0
        write(energy_final,string(energy[i],"\n"))
        filexyz.write(xyz_final,clusters[i])
        file_xyz_loc=open(string(folder,count,"_cluster_relax.xyz"),"w")
        filexyz.write(file_xyz_loc,clusters[i])
        close(file_xyz_loc)
        count += 1
    end
end
close(energy_final)
close(xyz_final)

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

# Counting the number of structure below a threshold of energy (eV)
energy_threshold=2
count=0
for i=1:nb_structure
    if abs(energy[i]-min_energy) < 2
        count += 1
    end
end
# From that point we care only about the structures below the threshold
nb_structure = count

# Building PIV-like vector for each cluster
nb_atoms=size(clusters[1].names)[1]
n_mo=Int(nb_atoms/3)
size_mo=Int(n_mo*(n_mo-1)/2)
n_s=Int(nb_atoms/3*2)
size_s=Int(n_s*(n_s-1)/2)
size_mos=n_s*n_mo
size_total=Int(nb_atoms*(nb_atoms-1)/2)
matrix=zeros(nb_structure, size_total )
for cluster=1:nb_structure
    #================================================#
    # compute Mo-Mo part
    for i=1:n_mo-1
        for j=i+1:n_mo
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[j,k])
            end
            matrix[cluster,(i-1)*n_mo+j-1]=sqrt(dist)
        end
    end
    for i=1:size_mo-1
        for j=i+1:size_mo
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    # Compute Mo-S part
    for i=1:n_mo
        for j=1:n_s
            dist=0
            for k=1:3
                dist += (clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+1+(i-1)*n_s+(j-1)]=sqrt(dist)
        end
    end
    for i=size_mo+1:size_mo+size_mos
        for j=i+1:size_mo+size_mos+1
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
    #Compute S-S part
    test=1
    for i=1:n_s-1
        for j=i+1:n_s
            dist=0
            for k=1:3
                dist=(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])*(clusters[cluster].positions[i,k]-clusters[cluster].positions[n_mo+j,k])
            end
            matrix[cluster,size_mo+size_mos+test]=sqrt(dist)
            test+=1
        end
    end
    for i=size_mo+size_mos+1:size_total-1
        for j=i+1:size_total
            if matrix[cluster,i] < matrix[cluster,j]
                dist=matrix[cluster,i]
                matrix[cluster,i]=matrix[cluster,j]
                matrix[cluster,j]=dist
            end
        end
    end
    #===============================================#
end

# Striking structure too close to those already existant
cut_off_distance=1.5 # Cut_off for topological distance, in Angstrom
strike=zeros(nb_structure)
for i=1:nb_structure-1
    for j=i+1:nb_structure
        dist=0
        for k=1:size_total
            dist += sqrt((matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k])*(matrix[i,k]-matrix[j,k]))
        end
        if dist/size_total < cut_off_distance
            strike[j]=1
        end
    end
end

# Writting remaining structures
energy_final=open(string(folder,"energy_final.dat"),"w")
xyz_final=open(string(folder,"all_relax_cluster.xyz"),"w")
count=1
for i=1:nb_structure
    if strike[i] == 0
        write(energy_final,string(energy[i],"\n"))
        filexyz.write(xyz_final,clusters[i])
        file_xyz_loc=open(string(folder,count,"_cluster_relax.xyz"),"w")
        filexyz.write(file_xyz_loc,clusters[i])
        close(file_xyz_loc)
        count += 1
    end
end
close(energy_final)
close(xyz_final)

#===============================================================================#
