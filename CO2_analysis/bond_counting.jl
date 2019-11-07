GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Counts the number of bonds (and averages over time windows)

using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using cpmd
using markov

Temperatures=[3000]
Volumes=[8.82]

stride_sim=5
fs2ps=0.001
time_sim=0.5 # in fs
unit=time_sim*fs2ps*stride_sim# in ps

nbC=32
nbO=64
nb_atoms=nbC+nbO

for T in Temperature

    for V in Volumes

        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
        file_in=string(folder,"TRAJEC_wrapped.xyz")
        folder_out=string(folder,"Data/")

        if ! isfile(file_in)
            continue
        end

        # Read trajectory
        atoms = filexyz.read( file_in )
        cell=cell_mod.Cell_param( V, V, V )

        nb_steps=size(atoms)[1]

        cut_off = [1.75]

        for co in cut_off

            nb_bonds=zeros(nb_steps)

            out=open(string(folder_out,"bonds_count_",V,"_",T,"_",co,".dat"),"w")

            # Computing number of bonds at each timestep
            for step=1:nb_steps
                print("progress : ",step/nb_steps*100,"%\n")
                count=0
                for carbon=1:32
                    for oxygen=33:96
                        if cell_mod.distance(atoms[step],cell,carbon,oxygen) < co
                            count += 1
                        end
                    end
                end
                nb_bonds[step]=count
                write(out,string(step*unit," ",count,"\n"))
            end
            close(out)

            total_sim=nb_steps*unit
            time_windows=[5]

            # Average over time windows
            for timewindow in time_windows

                nb_window=Int(trunc(total_sim/timewindow)+1)

                hist1D=zeros(nb_window)

                for step=1:nb_steps
                    for i=1:nb_window
                        if step*unit > (i-0.5)*timewindow && step*unit < (i+0.5)*timewindow
                            hist1D[i]+= nb_bonds[step]
                        end
                    end
                end

                out_box=open(string("/home/moogmt/bonds_count_",V,"_",T,"_",co,"A_box-",timewindow,".dat"),"w")
                for i=1:nb_window
                    write(out_box,string(i*timewindow," ",hist1D[i]*unit/timewindow,"\n"))
                end
                close(out_box)

            end
        end

    end
end
