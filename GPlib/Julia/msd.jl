include("contactmatrix.jl")

func="PBE"
T=3000
V=8.82

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
file="TRAJEC.xyz"

traj = filexyz.read(string(folder,file),stride)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

stop_100=20000

count=0

x0_std_block=traj[0].positions

file_std=open(string(folder,"MSD.dat"),"w")
file_block=open(string(folder,"MSD_block.dat"),"w")

for step=1:stop_100
    print("Progress: ",step/stop_100*100,"% \n")
    if step > nb_steps
        break
    end
    # WRAP
    msd_local_std=0
    msd_local_block=0
    for atom=1:nb_atoms
        for i=1:3
            msd_local_std += (traj[step].positions[atom,i]-traj[0].positions[atom,i])*(traj[step].positions[atom,i]-traj[0].positions[atom,i])
            msd_local_block += (traj[step].positions[atom,i]-x0_std_block[atom,i])*(traj[step].positions[atom,i]-x0_std_block[atom,i])
        end
    end
    msd_local /= nb_atoms
    msd_local_block /= nb_atoms
    count += 1
    if count == 1000
        x0_std_block = traj[step].positions
        count=0
    end
    msd_std +=  msd_local_std
    msd_block += msd_local_block
    write(file_std,string(step," ",msd_std,"\n"))
    write(file_block,string(step," ",msd_block,"\n"))
end

close(file_std)
close(file_block)
