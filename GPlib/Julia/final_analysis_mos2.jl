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
for n=2:4
    n_mo  = n
    n_s = 2*n
    file_all_data=open(string(folder_base,"Mo",n_mo,"S",n_s,"_unsorted_allinfo.dat"),"w")
    # 1.2.1 - Loop over all clusters
    file_cluster_all=open(string(folder_base,"cluster_centered.xyz"),"w")
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
                        else
                            print("n_mo: ",n_mo," n_s: ",n_s," elements:",elements,"\n")
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
            write(file_all_data,string(cluster," ",be_cluster," ",mag_abs," ",mag_tot,"\n"))
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
end

# Indicate possible redundancy for manual check
# Compute Absolute and Total Magnetization
# Compute HL Gap
