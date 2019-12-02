include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Time values
unit=0.005

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
                count_2 = stride
                count_3=1
                for step=1:nb_steps-1
                    print("Reading bond function - V: ",V," T: ",T,"K Progress: ",step/(nb_steps-1)*100,"%\n")
                    line=split(readline(file_bookeep))
                    if count_2 == stride
                        count_2=0
                        count=2
                        for carbon=1:nbC
                            for oxygen=1:nbO
                                bond_matrix[carbon,oxygen,count_3]=parse(Float64,line[count])
                                count += 1
                            end
                        end
                        #print("count: ",count_3,"\n")
                        count_3 += 1
                    end
                    count_2 += 1
                end
                close(file_bookeep)
            end
            #==================================================================#

            # clear memory
            atoms=[]

            # Carbon-Oxygen bonding
            #==================================================================#
            file_lifebonds=open(string(folder_out,"lifes_bonds.dat"),"w")
            for carbon=1:nbC
                for oxygen=1:nbO
                    # Not caring about stuff that do not move
                    check=sum(bond_matrix[carbon,oxygen,:])
                    if check == 0
                        continue
                    end
                    if check == nb_steps_analysis
                        write(file_lifebonds,string(carbon," ",oxygen," ",1," ",nb_steps_analysis," ",nb_steps_analysis," ",0,"\n"))
                        continue
                    end
                    step=1
                    while step < nb_steps_analysis
                        step_bonds=sum(bond_matrix[carbon,oxygen,step:nb_steps_analysis])
                        if step_bonds == 0
                            step = nb_steps_analysis
                        elseif bond_matrix[carbon,oxygen,step] == 0
                            step += 1
                        else
                            #print("V: ",V," T: ",T,"K - Counting lifetimes - C:",carbon," O:",oxygen," step: ",step,"\n")
                            life=1
                            write(file_lifebonds,string(carbon," ",oxygen," ",step," "))
                            #
                            stop_step=step
                            if step_bonds == nb_steps_analysis-step
                                write(file_lifebonds,string(nb_steps_analysis," ",nb_steps_analysis-step," ",0,"\n"))
                                step = nb_steps_analysis
                            else
                                for step2=step+1:nb_steps_analysis
                                    if bond_matrix[carbon,oxygen,step2] == 0
                                        stop_step=step2
                                    else
                                        stop_step += 1
                                        life += 1
                                    end
                                end
                                write(file_lifebonds,string(stop_step," ",life," "))
                                if stop_step < nb_steps_analysis
                                    write(file_lifebonds,string(1,"\n"))
                                else
                                    write(file_lifebonds,string(0,"\n"))
                                    break
                                end
                                step=stop_step
                            end
                        end
                    end
                end
            end
            close(file_lifebonds)
            #==================================================================#

        end
    end
end
