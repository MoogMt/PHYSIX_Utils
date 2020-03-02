# Compute the number of occurences of trimers
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

nbC=32
nbO=64

cut_off=1.75

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_base=string("/media/mathieu/Elements/CO2/")

for T in Temperatures
    for V in Volumes
        folder_in=string(folder_base,V,"/",T,"K/")
        folder_out=string( folder_in, "/Data/Trimer/" )
        file_traj = string( folder_in, "TRAJEC_fdb_wrapped.xyz" )
        if ! isfile( file_traj )
            continue
        end
        if !isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end
        cell=cell_mod.Cell_param(V,V,V)
        if ! isfile( file_traj )
            continue
        end
        traj = filexyz.readFileAtomList( file_traj )
        nb_steps = size( traj )[1]
        nb_atoms = size( traj[1].names)[1]
        handle_out = open( string( folder_out, "trimers_time.dat" ), "w" )
        for step=1:nb_steps
            print("Volume :",V," Temperature: ",T,"K Progress: ",round(step/nb_steps*100,digits=3),"% \n")
            cm = contact_matrix.buildMatrix( traj[step] , cell, cut_off )
            for carbon1=1:nbC-1
                for carbon2=carbon1+1:nbC
                    for oxygen1=1:nbO
                        if cm[ carbon1, nbC+oxygen1 ] > 0  && cm[ carbon2, nbC+oxygen1 ] > 0
                            for carbon3=carbon2+1:nbC
                                for oxygen2=1:nbO
                                    if cm[ carbon2, nbC+oxygen2 ] > 0  && cm[ carbon3, nbC+oxygen2 ] > 0
                                        for oxygen3=1:nbO
                                            if cm[ carbon1, nbC+oxygen3 ] > 0 && cm[ carbon3, nbC+oxygen3 ] > 0
                                                write( handle_out, string( step, " " ) )
                                                write( handle_out, string( carbon1, " " ) )
                                                write( handle_out, string( carbon2, " " ) )
                                                write( handle_out, string( carbon3, " " ) )
                                                write( handle_out, string( oxygen1+nbC, " " ) )
                                                write( handle_out, string( oxygen2+nbC, " " ) )
                                                write( handle_out, string( oxygen3+nbC, "\n") )
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        close( handle_out )
    end
end
