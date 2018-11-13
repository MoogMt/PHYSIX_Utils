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
    for cluster=1:n_max_cluster
        folder_local=string(folder_base,"Mo",n_mo,"S",n_s,"/Final/",cluster,"_cluster_relax/")
        if isfile(string(folder_local,"output"))
            print("Working on: Mo-",n_mo," S-",n_s," - cluster:",cluster,"\n")
            # 1.2.1.1 - Get Absolute Magnetization
            mag_hand=open(string(folder_local,"abs_mag_final"))
            lines=readlines(mag_hand)
            close(mag_hand)
            mag_abs=parse(Float64,lines[size(lines)[1]])
            # 1.2.1.2 - Get Absolute Magnetization
            # mag_hand_2=open(string(folder_local,"total_mag_final"))
            # lines=readlines(mag_hand_2)
            # close(mag_hand_2)
            # mag_total=parse(Float64,lines[size(lines)[1]])
            # 1.2.1.3 - Get Energy
            energy_hand = open(string(folder_local,"energy_final"))
            lines=readlines(energy_hand)
            close(energy_hand)
            e_cluster=parse(Float64,lines[size(lines)[1]])*Ry2eV
            # 1.2.1.4 - Compute Binding Energy
            nb_atoms=n_mo+n_s
            be_cluster = - ( e_cluster - n_mo*E_mo - n_s*E_s )/nb_atoms
            # 1.2.1.5 - Write Data
            write(file_all_data,string(cluster," ",be_cluster," ",mag_abs,"\n"))
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
            # priting cluster
            file_cluster=open(string(folder_local,"cluster_centered.xyz"),"w")
            write(file_cluster,string(nb_atoms,"\n"))
            write(file_cluster,string("Cluster: ",cluster,"\n"))
            for Mo=1:n_mo
                write(file_cluster,string("Mo "))
                for i=1:3
                    write(file_cluster,string(clusters[1].positions[Mo,i]," "))
                end
                write(file_cluster,string("\n"))
            end
            for S=1:n_s
                write(file_cluster,string("S "))
                for i=1:3
                    write(file_cluster,string(clusters[1].positions[n_mo+S,i]," "))
                end
                write(file_cluster,string("\n"))
            end
            close(file_cluster)
        end
    end
    close(file_all_data)
end

# Center each structure
# Indicate possible redundancy for manual check
# Compute Absolute and Total Magnetization
# Compute HL Gap
