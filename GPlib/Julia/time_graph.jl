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

# min number in which a bond has to be active for the autocorrelation to be done
min_steps=20

Volumes=[9.8]
Temperatures=[3000]

# for cut_off in Cut_Off
#     for T in Temperatures
#         for V in Volumes

            V=9.8
            T=3000
            cut_off=1.75

            folder_in=string(folder_base,V,"/",T,"K/")
            file=string(folder_in,"TRAJEC_wrapped.xyz")

            if ! isfile(string(folder_in,"TRAJEC_wrapped.xyz"))
                #continue
            end

            folder_out=string(folder_in,"Data/")

            atoms = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            nb_atoms=size(atoms[1].names)[1]
            nb_steps=size(atoms)[1]

            stride=5
            nb_steps_analysis=Int(trunc(nb_steps/stride))
            max_steps_autocor=Int(trunc(nb_steps/stride*frac))

            #==========================#
            # Bond matrix computation
            #==================================================================#
            bond_matrix=zeros(nbC,nbO,nb_steps_analysis)
            if ! isfile(string(folder_out,"bonds_book.dat"))
                file_bookeep=open(string(folder_out,"bonds_book.dat"),"w")
                for step=1:stride:nb_steps
                    print("Computing bond function - V: ",V," T: ",T,"K Progress: ",step/nb_steps*100,"%\n")
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

            #======================#
            # Bond Autocorrelation
            #==================================================================#
            time_corr_bonds=zeros(nbC,nbO,max_steps_autocor)
            for carbon=1:nbC
                for oxygen=1:nbO
                    # Determination of the number of steps for the autocorrelation
                    steps_corr=Int(sum(bond_matrix[carbon,oxygen,:]))
                    if steps_corr < min_steps
                        # If there is no bond except for flickering, we just don't care
                        continue
                    elseif steps_corr > nb_steps - min_steps
                        # If the bond is always there, no sens in doing autocorrelation
                        time_corr_bonds[carbon,oxygen,:] += 1
                        continue
                    end
                    # Boundary on the number of steps of autocorrelation to avoid errors
                    steps_corr = max_steps_autocor
                    # Autocorrelation
                    for i=1:steps_corr
                        print("V: ",V," T:",T,"K Autocorrelation - C:",carbon," O:",oxygen," Progress: ",i/(steps_corr-1)*100,"%\n")
                        for j=1:nb_steps_analysis-i+1
                            time_corr_bonds[carbon,oxygen,i] += bond_matrix[carbon,oxygen,j]*bond_matrix[carbon,oxygen,j+i-1]
                        end
                        time_corr_bonds[carbon,oxygen,i] /= (nb_steps_analysis-i+1)
                        time_corr_bonds[carbon,oxygen,i] /= time_corr_bonds[carbon,oxygen,1]
                    end
                end
            end
            #==================================================================#

            # Writting data to file
            #==================================================================#
            file_all=open(string(folder_out,"bond_correlation_individual.dat"),"w")
            for i=1:max_steps_autocor
                print("V:",V," T: ",T,"K Writting to file - Progress: ",i/(max_steps_autocor)*100,"%\n")
                write(file_all,string(i*unit*stride," "))
                for carbon=1:nbC
                    for oxygen=1:nbO
                        write(file_all,string(time_corr_bonds[carbon,oxygen,i]," "))
                    end
                end
                write(file_all,string("\n"))
            end
            close(file_all)
            #==================================================================#

            # Making histogram
            #==================================================================#
            time_int = 1 # ps
            nb_box_time = Int(trunc(nb_steps*unit/time_int))
            nb_box_corr=50
            inter_corr=1/nb_box_corr
            hist2d=zeros(nb_box_time,nb_box_corr)
            nb_steps=Int(trunc(nb_steps_analysis/nb_box_time)*nb_box_time)
            time_inter=Int(nb_steps/nb_box_time)
            for carbon=1:nbC
                for oxygen=1:nbO
                    print("Histogram Printing - C:",carbon," O:",oxygen,"\n")
                    for box=1:time_inter:nb_steps
                        for subbox=1:time_inter-1
                            for step=1:max_steps_autocor
                                for box2=1:nb_box_corr
                                    if time_corr_bonds < box2*inter_corr
                                        hist2d[box,box2] += 1
                                        break
                                    end
                                end
                            end
                        end
                    end
                end
            end
            # Printint histogram
            file_hist=open(string(folder_out,"bond_correlation_density.dat"),"w")
            for box_time=1:nb_box_time
                for box_corr=1:nb_box_corr
                    write(file_hist,string(box_time*unit*stride," ",box_corr*inter_corr," ",hist2d[box_time,box_corr],"\n"))
                end
                write(file_hist,string("\n"))
            end
            close(file_hist)
            #==================================================================#

#         end
#     end
# end
