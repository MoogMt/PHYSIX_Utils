include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]
Cut_Off=[1.75]

# Time values
unit=0.005
frac=0.5

# Number of atoms
nbC=32
nbO=nbC*2


Volumes=[8.82,9.8]
Temperatures=[2000,3000]

for cut_off in Cut_Off
    for T in Temperatures
        for V in Volumes


            folder_in=string(folder_base,V,"/",T,"K/")
            file=string(folder_in,"TRAJEC_wrapped.xyz")

            if ! isfile(folder_in,"TRAJEC_wrapped.xyz")
                continue
            end

            folder_out=string(folder_in,"Data/")

            atoms = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            nb_atoms=size(atoms[1].names)[1]
            nb_steps=size(atoms)[1]
            steps=Int(trunc(nb_steps*frac))

            bond_matrix=zeros(nbC,nbO,nb_steps)
            if ! isfile(string(folder_out,"bonds_book.dat"))
                file_bookeep=open(string(folder_out,"bonds_book.dat"),"w")
                for step=1:nb_steps
                    print("V: ",V," T: ",T,"K Progress: ",step/nb_steps*100,"%\n")
                    write(file_bookeep,string(step," "))
                    for carbon=1:nbC
                        for oxygen=1:nbO
                            out=0
                            if cell_mod.distance(atoms[step],cell,carbon,nbC+oxygen) < cut_off
                                out += 1
                            end
                            bond_matrix[carbon,oxygen,step]=out
                            write(file_bookeep,string(out," "))
                        end
                    end
                    write(file_bookeep,string("\n"))
                end
                close(file_bookeep)
            else
                file_bookeep=open(string(folder_out,"bonds_book.dat"))
                for step=1:nb_steps
                    line=split(readline(file_bookeep))
                    count=2
                    for carbon=1:nbC
                        for oxygen=1:nbO
                            bond_matrix[carbon,oxygen,step]=parse(Float64,line[count])
                            count += 1
                        end
                    end
                end
                close(file_bookeep)
            end

            time_corr_bonds=zeros(nbC,nbO,steps)
            for carbon=1:nbC
                for oxygen=1:nbO
                    # If there is no bond except for flickering, we just don't care
                    if sum(bond_matrix[carbon,oxygen,:]) < 20
                        continue
                    end
                    # Autocorrelation
                    for i=1:steps
                        print("V: ",V," T:",T,"K Autocorrelation - C:",carbon," O:",oxygen," Progress: ",i/(steps-1)*100,"%\n")
                        for j=1:nb_steps-i+1
                            time_corr_bonds[carbon,oxygen,i] += bond_matrix[carbon,oxygen,j]*bond_matrix[carbon,oxygen,j+i-1]
                        end
                        time_corr_bonds[carbon,oxygen,i] /= (nb_steps-i+1)
                    end
                end
            end

            file_all=open(string(folder_out,"bond_correlation_individual.dat"),"w")
            for i=1:steps
                print("V:",V," T: ",T,"K Writting to file - Progress: ",i/(steps)*100,"%\n")
                write(file_all,string(i*unit," "))
                for carbon=1:nbC
                    for oxygen=1:nbO
                        write(file_all,string(time_corr_bonds[carbon,oxygen,i]," "))
                    end
                end
                write(file_all,string("\n"))
            end
            close(file_all)

        end
    end
end
