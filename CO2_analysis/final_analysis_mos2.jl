include("contactmatrix.jl")

Bohr2Ang=0.529177
Ry2eV=13.6056980659

folder_base="/media/moogmt/Stock/MoS2/Relax/PBE/"

n_max_cluster=200

# 1 - Compute BE of each structure

# 1.1 - Getting Mo and S energy and conversion from Ry to eV

file_mo=string("/media/moogmt/Stock/MoS2/Relax/PBE/Mo/30B/energy_final")
fileh_mo=open(file_mo)
e_mo=parse(Float64,readlines(fileh_mo)[1])*Ry2eV
close(fileh_mo)

file_s=string("/media/moogmt/Stock/MoS2/Relax/PBE/S/30B/energy_final")
fileh_s=open(file_s)
e_s=parse(Float64,readlines(fileh_s)[1])*Ry2eV
close(fileh_s)

a=30*Bohr2Ang

# 1.2 Loop over all sizes
file_energy=open(string(folder_base,"BE-n.dat"),"w")
for n=2:4
    n_mo  = n
    n_s = 2*n
    nb_atoms=n_mo+n_s
    file_all_data=open(string(folder_base,"Mo",n_mo,"S",n_s,"_unsorted_allinfo.dat"),"w")
    # 1.2.1 - Loop over all clusters
    file_cluster_all=open(string(folder_base,"Mo",n_mo,"S",n_s,"cluster_centered.xyz"),"w")
    for cluster=1:n_max_cluster
        folder_local=string(folder_base,"Mo",n_mo,"S",n_s,"/Final/",cluster,"_cluster_relax/")
        if isfile(string(folder_local,"output"))
            print("Working on: Mo-",n_mo," S-",n_s," - cluster:",cluster,"\n")
            # 1.2.1.1 - Get Electronic stuff
            file_output=open(string(folder_local,"output"))
            lines=readlines(file_output)
            close(file_output)
            mag_tot=0
            mag_abs=0
            hl_gap=0
            fermi_energy=0
            for line_nb=1:size(lines)[1]
                if ! isempty(lines[line_nb])
                    elements=split(lines[line_nb])
                    # HL Gap
                    if size(elements)[1] < 3
                        continue
                    end
                    if elements[2] == "SPIN" && elements[3] == "UP"
                        bands=[]
                        start=line_nb+5
                        while ! isempty( lines[start] )
                            elements2=split(lines[start])
                            for i=1:size(elements2)[1]
                                push!(bands,parse(Float64,elements2[i]))
                            end
                            start += 1
                        end
                        start += 6
                        while ! isempty( lines[start] )
                            elements2=split(lines[start])
                            for i=1:size(elements2)[1]
                                push!(bands,parse(Float64,elements2[i]))
                            end
                            start += 1
                        end
                        fermi_energy=parse(Float64,split(lines[start+1])[5])
                    elseif elements[2] == "magnetization"
                        if elements[1] == "total"
                            mag_tot=parse(Float64,elements[4])
                        elseif  elements[5] == "Bohr"
                            mag_abs=parse(Float64,elements[4])
                        end
                    end
                end
            end
            bands -= fermi_energy
            min=0
            for i=1:size(bands)[1]
                if bands[i] > 0
                    if min < bands[i]
                        min=bands[i]
                    end
                end
            end
            max=0
            for i=1:size(bands)[1]
                if bands[i] > 0
                    if max > bands[i]
                        max=bands[i]
                    end
                end
            end
            # Compute HL Gap
            hl_gap=min+max
            #-----------------------------------------------------------------
            # 1.2.1.3 - Get Energy
            energy_hand = open(string(folder_local,"energy_final"))
            lines=readlines(energy_hand)
            close(energy_hand)
            e_cluster=parse(Float64,lines[size(lines)[1]])*Ry2eV
            # 1.2.1.4 - Compute Binding Energy
            nb_atoms=n_mo+n_s
            be_cluster = - ( e_cluster - n_mo*e_mo - n_s*e_s )/nb_atoms
            #-----------------------------------------------------------------
            # 1.2.1.5 - Write Data
            write(file_all_data,string(cluster," ",be_cluster," ",mag_abs," ",mag_tot," ",hl_gap,"\n"))
            #-----------------------------------------------------------------
            # 1.2.1.6 - Centering cluster
            clusters=filexyz.readFastFile(string(folder_local,"cluster.xyz"))
            # Computing barycenter
            bc=zeros(3)
            for atom=1:nb_atoms
                for i=1:3
                    bc[i] += clusters[1].positions[atom,i]
                end
            end
            bc /= nb_atoms
            # Centering cluster
            center=ones(3)*a*0.5
            move_vector=bc-center
            for atom=1:nb_atoms
                clusters[1].positions[atom,:] -= move_vector
            end
            # Wrapping
            for atom=1:nb_atoms
                for i=1:3
                    clusters[1].positions[atom,i] = cell_mod.wrap(clusters[1].positions[atom,i],a)
                end
            end
            # Re-Centering cluster
            center=ones(3)*a*0.5
            move_vector=bc-center
            for atom=1:nb_atoms
                clusters[1].positions[atom,:] -= move_vector
            end
            # priting cluster
            file_cluster=open(string(folder_local,"cluster_centered.xyz"),"w")
            write(file_cluster,string(nb_atoms,"\n"))
            write(file_cluster_all,string(nb_atoms,"\n"))
            write(file_cluster,string("Cluster: ",cluster,"\n"))
            write(file_cluster_all,string("Cluster: ",cluster,"\n"))
            for Mo=1:n_mo
                write(file_cluster,string("Mo "))
                write(file_cluster_all,string("Mo "))
                for i=1:3
                    write(file_cluster,string(clusters[1].positions[Mo,i]," "))
                    write(file_cluster_all,string(clusters[1].positions[Mo,i]," "))
                end
                write(file_cluster,string("\n"))
                write(file_cluster_all,string("\n"))
            end
            for S=1:n_s
                write(file_cluster,string("S "))
                write(file_cluster_all,string("S "))
                for i=1:3
                    write(file_cluster,string(clusters[1].positions[n_mo+S,i]," "))
                    write(file_cluster_all,string(clusters[1].positions[n_mo+S,i]," "))
                end
                write(file_cluster,string("\n"))
                write(file_cluster_all,string("\n"))
            end
            close(file_cluster)
        end
    end
    close(file_all_data)
    close(file_cluster_all)

    # Sorting clusters
    file_all_data=open(string(folder_base,"Mo",n_mo,"S",n_s,"_unsorted_allinfo.dat"))
    lines_data=readlines(file_all_data)
    close(file_all_data)
    be=zeros(size(lines_data)[1])
    hl_gap=zeros(size(lines_data)[1])
    mag_abs=zeros(size(lines_data)[1])
    mag_tot=zeros(size(lines_data)[1])
    index=zeros(size(lines_data)[1])
    for i=1:size(lines_data)[1]
        elements=split(lines_data[i])
        be[i]=parse(Float64,elements[2])
        mag_abs[i]=parse(Float64,elements[3])
        mag_tot[i]=parse(Float64,elements[4])
        hl_gap[i]=parse(Float64,elements[5])
        index[i]=i
    end
    clusters=filexyz.readFastFile(string(folder_base,"Mo",n_mo,"S",n_s,"cluster_centered.xyz"))
    for i=1:size(clusters)[1]
        for j=1:size(clusters)[1]
            if be[i] > be[j]
                #--------------------------
                stock=be[i]
                be[i]=be[j]
                be[j]=stock
                #-------------------------
                stock=mag_abs[i]
                mag_abs[i]=mag_abs[j]
                mag_abs[j]=stock
                #------------------------
                stock=mag_tot[i]
                mag_tot[i]=mag_tot[j]
                mag_tot[j]=stock
                #-----------------------
                # Mag tot
                #--------------------------
                stock=hl_gap[i]
                hl_gap[i]=hl_gap[j]
                hl_gap[j]=stock
                #--------------------------
                stock=clusters[i]
                clusters[i]=clusters[j]
                clusters[j]=stock
                #--------------------------
                stock=index[i]
                index[i]=index[j]
                index[j]=stock
                #--------------------------
            end
        end
    end
    strike=zeros(size(lines_data)[1])
    for i=1:size(clusters)[1]-1
        if strike[i] == 1
            continue
        end
        for j=i+1:size(clusters)[1]
            if abs(be[i]-be[j]) < 0.001
                strike[j]=1
            end
        end
    end
    file_all_data=open(string(folder_base,"Mo",n_mo,"S",n_s,"_sorted_allinfo.dat"),"w")
    file_all_cluster=open(string(folder_base,"Mo",n_mo,"S",n_s,"_all_sorted_clusters.xyz"),"w")
    file_cluster=open(string(folder_base,"Mo",n_mo,"S",n_s,"_cluster",count,".xyz"),"w")
    count=1
    for i=1:size(clusters)[1]
        if strike[i] == 1
            print("STRIKE: ",i,"\n")
            continue
        end
        file_cluster=open(string(folder_base,"Mo",n_mo,"S",n_s,"_cluster",count,".xyz"),"w")
        write(file_all_data,string(count," ",be[i]," ",mag_abs[i]," ",mag_tot[i]," ",hl_gap[i]," ",index[i],"\n"))
        write(file_all_cluster,string(nb_atoms,"\n"))
        write(file_all_cluster,string("Cluster: ",index[i],"\n"))
        write(file_cluster,string(nb_atoms,"\n"))
        write(file_cluster,string("Cluster: ",count,"\n"))
        for Mo=1:n_mo
            write(file_all_cluster,string("Mo "))
            write(file_cluster,string("Mo "))
            for j=1:3
                write(file_all_cluster,string(clusters[i].positions[Mo,j]," "))
                write(file_cluster,string(clusters[i].positions[Mo,j]," "))
            end
            write(file_all_cluster,string("\n"))
            write(file_cluster,string("\n"))
        end
        for S=1:n_s
            write(file_all_cluster,string("S "))
            write(file_cluster,string("S "))
            for j=1:3
                write(file_all_cluster,string(clusters[i].positions[n_mo+S,j]," "))
                write(file_cluster,string(clusters[i].positions[n_mo+S,j]," "))
            end
            write(file_all_cluster,string("\n"))
            write(file_cluster,string("\n"))
        end
        count+=1
        close(file_cluster)
    end
    close(file_all_data)
    close(file_all_cluster)
    write(file_energy,string(n," ",be[1],"\n"))
end

close(file_energy)
