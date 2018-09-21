include("contactmatrix.jl")

func="PBE-MT"
Volumes=[8.82]
Temperatures=[2000,2250,2500,3000]
stride=1
unit=0.005
block_size=[500,1000]

nbC = 32
nbO = 64

for size_block in block_size

    for V in Volumes

        for T in Temperatures

            folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
            file="TRAJEC.xyz"

            if isfile( string(folder,file) )

                traj = filexyz.read( string(folder,file), stride )
                cell=cell_mod.Cell_param( V, V, V )

                nb_steps=size( traj )[1]
                nb_atoms=size( traj[1].names )[1]

                stop_100 = 20000

                count=0
                msd_stock   = 0
                msd_stock_C = 0
                msd_stock_O = 0

                bary_origin   = zeros(3)
                bary_origin_C = zeros(3)
                bary_origin_O = zeros(3)
                bary_step   = zeros(3)
                bary_step_C = zeros(3)
                bary_step_O = zeros(3)

                x0_std_block = traj[1].positions[:,:]

                # Barycenter computation
                #==============================================================#
                for atom=1:nb_atoms
                    for i=1:3
                        bary_origin[i] += traj[1].positions[atom,i]
                        if atom < 33
                            bary_origin_C[i] += traj[1].positions[atom,i]
                        else
                            bary_origin_O[i] += traj[1].positions[atom,i]
                        end
                    end
                end
                for i=1:3
                    bary_origin[i]   /= nb_atoms
                    bary_origin_C[i] /= nbC
                    bary_origin_O[i] /= nbO
                end
                #==============================================================#

                #==============================================================#
                file_std=open(string(folder,"MSD.dat"),"w")
                file_std_C=open(string(folder,"MSD_C.dat"),"w")
                file_std_O=open(string(folder,"MSD_O.dat"),"w")
                file_block=open(string(folder,"MSD_block-",size_block,".dat"),"w")
                file_block_C=open(string(folder,"MSD_C_block-",size_block,".dat"),"w")
                file_block_O=open(string(folder,"MSD_O_block-",size_block,".dat"),"w")
                #==============================================================#

                # Step
                #==============================================================#
                for step=1:stop_100
                    #==============================================================#
                    print("V = ",V," T = ",T," Progress: ",step/stop_100*100,"% \n")
                    #==============================================================#

                    #==============================================================#
                    if step > nb_steps
                        break
                    end
                    #==============================================================#

                    #==============================================================#
                    bary_step   = zeros(3)
                    bary_step_C = zeros(3)
                    bary_step_O = zeros(3)
                    for atom=1:nb_atoms
                        for i=1:3
                            bary_step[i] += traj[step].positions[atom,i]
                            if atom < 33
                                bary_step_C[i] += traj[step].positions[atom,i]
                            else
                                bary_step_O[i] += traj[step].positions[atom,i]
                            end
                        end
                    end
                    for i=1:3
                        bary_step[i]  /= nb_atoms
                        bary_step_C[i] /= nbC
                        bary_step_O[i] /= nbO
                    end
                    #==============================================================#

                    # WRAP
                    #==========================================================#
                    msd_local_std = 0
                    msd_local_block = 0
                    msd_local_std_C = 0
                    msd_local_std_O = 0
                    msd_local_block_C = 0
                    msd_local_block_O = 0
                    for atom=1:nb_atoms
                        for i=1:3
                            msd_local_std += ( ( traj[step].positions[atom,i] - bary_step[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )
                            msd_local_block += ( ( traj[step].positions[atom,i] - bary_step[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )
                            if atom < 33
                                msd_local_std_C += ( ( traj[step].positions[atom,i] - bary_step_C[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step_C[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )
                                msd_local_block_C += ( ( traj[step].positions[atom,i] - bary_step_C[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step_C[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )
                            else
                                msd_local_std_O += ( ( traj[step].positions[atom,i] - bary_step_O[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step_O[i] ) - ( traj[1].positions[atom,i] - bary_origin[i] ) )
                                msd_local_block_O += ( ( traj[step].positions[atom,i] - bary_step_O[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )*( ( traj[step].positions[atom,i] - bary_step_O[i] ) - ( x0_std_block[atom,i] - bary_origin[i] ) )
                            end
                        end
                    end
                    #==========================================================#

                    #==========================================================#
                    write( file_std,   string( step*unit," ", msd_local_std/nb_atoms,"\n"   ) )
                    write( file_std_C, string( step*unit," ", msd_local_std_C/nb_atoms,"\n" ) )
                    write( file_std_O, string( step*unit," ", msd_local_std_O/nb_atoms,"\n" ) )
                    #==========================================================#

                    #==========================================================#
                    msd_local_block = msd_local_block/nb_atoms+ msd_stock
                    msd_local_block_C = msd_local_block_C/nb_atoms + msd_stock_C
                    msd_local_block_O = msd_local_block_O/nb_atoms + msd_stock_O
                    write( file_block,   string( step*unit," ",msd_local_block,"\n"   ) )
                    write( file_block_C, string( step*unit," ",msd_local_block_C,"\n" ) )
                    write( file_block_O, string( step*unit," ",msd_local_block_O,"\n" ) )
                    #==========================================================#

                    #==========================================================#
                    if count == block_size
                        x0_std_block = traj[step].positions
                        msd_stock = msd_local_block
                        msd_stock_C = msd_local_block_C
                        msd_stock_O = msd_local_block_O
                        count=0
                    end
                    count += 1
                    #==========================================================#

                end
                #==============================================================#

                #==========================================================#
                close( file_std )
                close( file_std_C )
                close( file_std_O )
                close( file_block )
                close( file_block_C )
                close( file_block_O )
                #==========================================================#

            end

        end

    end

end
