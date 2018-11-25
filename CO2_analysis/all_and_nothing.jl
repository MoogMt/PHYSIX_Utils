include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

for cut_off in Cut_Off
    for T in Temperatures
        for V in Volumes

            folder_in=string(folder_base,V,"/",T,"K/")
            file=string(folder_in,"TRAJEC_wrapped.xyz")

            if ! isfile(string(folder_in,"TRAJEC_wrapped.xyz"))
                continue
            end

            folder_out=string(folder_in,"Data/")

            atoms = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            nb_atoms=size(atoms[1].names)[1]
            nb_steps=size(atoms)[1]

            stride=5
            nb_steps_analysis=Int(trunc(nb_steps/stride))
            max_steps_autocor=Int(trunc(nb_steps/stride*frac))

            restart_bonds_book = false

            #==========================#
            # Bond matrix computation
            #==================================================================#
            bond_matrix=zeros(nbC,nbO,nb_steps_analysis)
            if ( ! isfile(string(folder_out,"bonds_book.dat")) )|| restart_bonds_book == true
                file_bookeep=open(string(folder_out,"bonds_book.dat"),"w")
                count=0
                step2=1
                for step=1:nb_steps
                    print("Computing bond function - V: ",V," T: ",T,"K Progress: ",step/nb_steps*100,"%\n")
                    write(file_bookeep,string(step," "))
                    check=false
                    for carbon=1:nbC
                        for oxygen=1:nbO
                            out=0
                            if cell_mod.distance(atoms[step],cell,carbon,nbC+oxygen) < cut_off
                                out += 1
                            end
                            if count == stride
                                bond_matrix[carbon,oxygen,step2]=out
                                check=true
                            end
                            write(file_bookeep,string(out," "))
                        end
                    end
                    if check
                        step2 += 1
                        count=0
                    end
                    count+=1
                    write(file_bookeep,string("\n"))
                end
                close(file_bookeep)
            else
                file_bookeep=open(string(folder_out,"bonds_book.dat"))
                for step=1:stride:nb_steps_analysis
                    print("Reading bond function - V: ",V," T: ",T,"K Progress: ",step/nb_steps_analysis*100,"%\n")
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
            #==================================================================#

            # clear memory
            atoms=[]

            file_all=open(string(folder_out,"bonds_all_book.dat"),"w")
            file_all_book=open(string(folder_out,"bonds_all_book.dat"),"w")
            file_nothing=open(string(folder_out,"bonds_nothing.dat"),"w")
            file_nothing_book=open(string(folder_out,"bonds_nothing_book.dat"),"w")
            count_all=0
            count_nothing=0
            for carbon=1:nbC
                for oxygen=1:nbO
                    print(" Counting - C: ",C," O:",O,"\n")
                    sum_of_all_bonds=sum(bond_matrix[carbon,oxygen,:])
                    if sum_of_all_bonds == 0
                        count_nothing += 1
                        write(file_nothing_book,string(carbon," ",oxygen,"\n"))
                    elseif sum_of_all_bonds == nb_steps_analysis
                        count_all += 1
                        write(file_all_book,string(carbon," ",oxygen,"\n"))
                    end
                end
            end
            write(file_all,string(count_all/nbC*nbO,"\n"))
            write(file_nothing,string(count_nothing/(nbC*nbO),"\n"))
            close(file_all)
            close(file_nothing)

        end
    end
end
