include("contactmatrix.jl")

func="PBE-MT"
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]


unit=0.005

nb_box=1

nbC = 32
nbO = 64

nb_box=25
box_time=200

for V in Volumes

    for T in Temperatures

        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
        file="TRAJEC.xyz"

        if isfile( string(folder,file) )
            traj = filexyz.readFastFile( string(folder,file))
            cell=cell_mod.Cell_param( V, V, V )
            nb_steps=size( traj )[1]
            nb_atoms=size( traj[1].names )[1]
            folder_out=string(folder,"Data/")
            file_std=open(string(folder_out,"MSD_O.dat"),"w")
            # Determining max and computing MSD
            msd_local_std = zeros(nb_steps,nb_atoms)
            max=0
            min=0
            for step=1:nb_steps
                msd=0
                print("Computing MSD - V = ",V," T = ",T," Progress: ",step/nb_steps*100,"% \n")
                write( file_std, string(step*unit," ",))
                for atom=1:nbO
                    for i=1:3
                        msd_local_std[step,atom] +=(traj[step].positions[nbC+atom,i]-traj[1].positions[nbC+atom,i])*(traj[step].positions[nbC+atom,i]-traj[1].positions[nbC+atom,i])
                    end
                    if msd_local_std[step,atom] > max
                        max=msd_local_std[step,atom]
                    end
                    write( file_std, string( msd_local_std[step,atom]," " ) )
                end
                write(file_std,string("\n"))
            end
            close( file_std )
            # Making histogram

            interval=(max-min)/nb_box
            file_hist=open(string(folder_out,"MSD_O_histogram.dat"),"w")
            count_step=1
            box_hist=zeros(nb_box)
            for step=1:nb_steps
                print("Histogram - V = ",V," T = ",T," Progress: ",step/nb_steps*100,"% \n")
                count=0
                for atom=1:nb_atoms
                    for box=1:nb_box
                        if msd_local_std[step,atom] > (box-1)*interval && msd_local_std[step,atom] < box*interval
                            box_hist[box] +=1
                            count +=1
                            break
                        end
                    end
                end
                box_hist /= count
                for box=1:nb_box
                    write(file_hist,string(count_step*unit*box_time," ",box*interval," ",box_hist[box],"\n"))
                end
                write(file_hist,string("\n"))
                if ( step % box_time ) == 0
                    box_hist=zeros(nb_box)
                    count_step+=1
                end
            end
            close(file_hist)
        end
    end
end
