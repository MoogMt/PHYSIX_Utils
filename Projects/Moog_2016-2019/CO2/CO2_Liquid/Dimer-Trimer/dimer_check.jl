# Loading file
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

nbC=32
nbO=64

folder_base = string( "/media/mathieu/Elements/CO2/" )

cut_off=1.75

Volumes      = [ 10.0, 9.8, 9.5, 9.4, 9.375, 9.35]#, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6 ]
Temperatures = [ 1750, 2000, 2500, 3000 ]

for T in Temperatures
    for V in Volumes
        folder_target = string( folder_base, V, "/", T, "K/" )
        file_traj = string( folder_target, "TRAJEC_fdb_wrapped.xyz" )
        if ! isfile( file_traj )
            continue
        end
        folder_out = string( folder_target, "Data/Dimers/" )
        if ! isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end
        traj = filexyz.readFileAtomList( file_traj )
        cell = cell_mod.Cell_param( V, V, V )
        nb_steps = size( traj )[1]
        nb_atoms = size( traj[1].names )[1]
        dimer_out = open( string( folder_out, "dimer-", cut_off, ".dat" ), "w" )
        for step=1:nb_steps
            print( "Volume :", V, " Temperature: ", T, "K Progress: ", round(step/nb_steps*100, digits=3 ), "% \n" )
            cm = contact_matrix.buildMatrix( traj[step], cell, cut_off )
            used = zeros( nb_atoms )
            for carbon1=1:nbC-1
                for carbon2=carbon1+1:nbC
                    for oxygen1=1:nbO
                        if cm[ carbon1, nbC+oxygen1 ] > 0  && cm[ carbon2, nbC+oxygen1 ] > 0
                            for oxygen2=oxygen1+1:nbO
                                if cm[ carbon1, oxygen2 ] > 0  && cm[ carbon2, nbC+oxygen2 ] > 0
                                    write( dimer_out, string( step, " ", carbon1, " ", carbon2," ", nbC+oxygen1, " ", nbC+oxygen2, "\n" ) )
                                end
                            end
                        end
                    end
                end
            end
        end
        close( dimer_out )
    end
end
