GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Computes the distances of the four nearest oxygen atoms to carbon atoms,
# Along with their angles

# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz

volumes=[8.82]
temperatures=[3000]

max_neigh=4

ps2fs=0.001
timestep=0.5
stride = 1
unit=ps2fs*timestep*stride

nbC=32
nbO=64
nb_atoms=nbC+nbO

folder_base=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

for V in volumes
    for T in temperatures

        folder=string(folder_base,V,"/",T,"K/")
        file=string(folder,"TRAJEC_wrapped.xyz")

        if ! isfile( string(file) )
            continue
        end

        # Reading trajectory
        traj = filexyz.read(file,stride,start_step)
        cell=cell_mod.Cell_param(V,V,V)

        nb_steps=size(traj)[1]

        # Folder for output data
        folder_out=string(folder_base,"Data/")

        fileC=open(string(folder_out,"distanglesC-",V,"-",T,"K.dat"),"w")
        # Loop over steps
        for step=1:nb_steps
            print("Progres: ",step/nb_steps*100,"%\n")
            # Loop over carbons
            for carbon=1:nbC
                # Keep distances of oxygen to carbon
                distances = zeros(nbO)
                for oxygen=1:nbO
                    distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
                end
                # get the index of the oxygen sorted by distances
                index=sortperm(distances)

                # write distances to file
                write(fileC, string(step*unit," ") )
                for i=1:max_neigh
                    write(fileC, string(distances[index[i]]," ") )
                end

                # Computes and write the angles between the N nearest oxygen
                for i=1:max_neigh-1
                    for j=i+1:max_neigh
                        a=cell_mod.distance(traj[step],cell,carbon,Int(index[i]+nbC))
                        b=cell_mod.distance(traj[step],cell,carbon,Int(index[j]+nbC))
                        c=cell_mod.distance(traj[step],cell,Int(index[i]+nbC),Int(index[j]+nbC))
                        write(fileC,string(acosd((a*a+b*b-c*c)/(2*a*b))," "))
                    end
                end

                write(fileC,string("\n"))
            end
        end
        close(fileC)

    end
end
