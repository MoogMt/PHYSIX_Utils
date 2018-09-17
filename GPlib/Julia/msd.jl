include("contactmatrix.jl")

func="PBE-MT"
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.8,10.0]
Temperatures=[2000,2250,2500,3000]

unit=0.005

for V in Volumes

    for T in Temperatures


        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
        file="TRAJEC.xyz"

        traj = filexyz.read(string(folder,file),stride)
        cell=cell_mod.Cell_param(V,V,V)

        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]

        stop_100=20000

        count=0

        x0_std_block=traj[1].positions

        file_std=open(string(folder,"MSD.dat"),"w")
        file_block=open(string(folder,"MSD_block.dat"),"w")

        for step=1:stop_100
            print("V = ",V," T = ",T," Progress: ",step/stop_100*100,"% \n")
            if step > nb_steps
                break
            end
            # WRAP
            msd_local_std = 0
            msd_local_block = 0
            for atom=1:nb_atoms
                for i=1:3
                    msd_local_std += (traj[step].positions[atom,i]-traj[1].positions[atom,i])*(traj[step].positions[atom,i]-traj[1].positions[atom,i])
                    msd_local_block += (traj[step].positions[atom,i]-x0_std_block[atom,i])*(traj[step].positions[atom,i]-x0_std_block[atom,i])
                end
            end
            count += 1
            if count == 1000
                x0_std_block = traj[step].positions
                count=0
            end
            msd_std +=  msd_local_std
            msd_block += msd_local_block
            write(file_std,string(step*unit," ",msd_local_std,"\n"))
            write(file_block,string(step*unit," ",msd_local_block,"\n"))
        end

        close(file_std)
        close(file_block)

    end

end
