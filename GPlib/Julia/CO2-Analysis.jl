include("contactmatrix.jl")

# Thermodynamical values
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

stride=1
unit=0.005

for cut_off in Cut_Off
    for volume in Volumes
        for temperature in Temperatures

            print("Processing : ",volume," ",temperature,"\n")

            folder=string(folder_base,volume,"/",temperature,"K/")

            file_traj = string(folder,"TRAJEC_wrapped.xyz")
            molecule_file_name = string(folder,"Data/molecules_allinfo_cutoff-",cut_off,".dat")
            molecule_size_file_name = string(folder,"Data/molecules_sizes_cutoff-",cut_off,".dat")

            if isfile( file_traj )

                # Reading file
                print("Reading xyz file","\n")
                traj = filexyz.readFastFile( file_traj )
                cell=cell_mod.Cell_param( volume, volume, volume )

                # Number of atoms
                nb_steps=size(traj)[1]
                nb_atoms=size(traj[1].names)[1]
                sizes=zeros(nb_atoms)

                # Computing molecule information
                #---------------------------------------------------------------
                print("Processing molecules: ",volume," ", temperature,"\n")
                file_molecule = open( molecule_file_name , "w" )
                for step=1:nb_steps

                    print("Processing molecules - ",volume," - ", temperature, " - ",step/nb_steps*100,"%\n")

                    # Creating bond matrix
                    #-----------------------------------------------------------
                    matrix_bonds=zeros( nb_atoms, nb_atoms )
                    for i=1:nb_atoms
                        for j=i+1:nb_atoms
                            if cell_mod.distance( traj[step], cell, i, j ) < cut_off
                                matrix_bonds[i,j] = 1
                                matrix_bonds[j,i] = 1
                            end
                        end
                    end
                    #-----------------------------------------------------------

                    # Counting the number of molecules\
                    #-----------------------------------------------------------
                    nb_mol=0
                    mol_index=zeros(nb_atoms)
                    for i=1:nb_atoms
                        if mol_index[i] == 0
                            nb_mol += 1
                            mol_index[i]=nb_mol
                            mol_index = graph_mod.searchGroupMember(matrix_bonds,mol_index,i,nb_mol)
                        end
                    end
                    #-----------------------------------------------------------

                    # Writting molecule information to disk
                    #-----------------------------------------------------------
                    for mol_index_loop=1:nb_mol
                        size=0
                        write(file_molecule,string(step," ",nb_mol," "))
                        write(file_molecule,string(mol_index_loop," "))
                        for atom=1:nb_atoms
                            if mol_index[atom] == mol_index_loop
                                size += 1
                                write(file_molecule,string(atom," "))
                            end
                        end
                        sizes[ size ] += 1
                        write(file_molecule,string(size," \n"))
                    end
                    #-----------------------------------------------------------

                end
                close(file_molecule)
                #-----------------------------------------------------------

                # Registring all size statistics
                #-----------------------------------------------------------
                file_size_molecule= open( molecule_size_file_name, "w" )
                for i=1:nb_atoms
                    write(file_size_molecule,string(i," ",sizes[i]," ",sizes[i]/nb_atoms*i,"\n"))
                end
                close(file_size_molecule)
                #-----------------------------------------------------------

                #-----------------------------------------------------------


            else
                print("File: ",file_traj," does not exists, moving on.\n")
            end

        end
    end
end
