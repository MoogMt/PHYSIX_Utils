include("contactmatrix.jl")

# Thermodynamical values
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.35,9.4,9.5,9.8]
Temperature=[2000,2250,2500,3000,3500]

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

stride=5
unit=0.0005

cut_off = 1.75

nb_steps=40000

for volume in Volumes
    for temperature in Temperatures

        folder=string(folder_base,current_volume,"/",temperature,"K/")
        file_traj=string(folder,"TRAJEC_wrapped.xyz")

        if isfile( file_traj )

            # Reading file
            traj = filexyz.readFastFile( file_traj )
            cell=cell_mod.Cell_param( volume, volume, volume )

            # Determining number of steps
            if nb_steps > size(atoms)[1]
                nb_steps=size(atoms)[1]
            end

            # Number of atoms
            nb_atoms=size(atoms[1].names)[1]

            # Computing molecule information
            #---------------------------------------------------------------
            file_molecule = open( string(folder,"molecules_allinfo_cutoff-",cut_off,".dat"), "w" )
            for step=1:nb_steps

                # Creating bond matrix
                #-----------------------------------------------------------
                matrix_bonds=zeros(nb_atoms,nb_atoms)
                for i=1:nb_atoms
                    for j=i+1:nb_atoms
                        if cell_mod.distance(atoms[step],cell,i,j) < cut_off
                            matrix_bonds[i,j]=1
                            matrix_bonds[j,i]=1
                        end
                    end
                end
                #-----------------------------------------------------------

            end
            close(file_molecule)
        else
            print("File: ",file_traj," does not exists, moving on.\n")
        end

    end
end
