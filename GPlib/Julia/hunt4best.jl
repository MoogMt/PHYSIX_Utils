include("contactmatrix.jl")

func="PBE-MT"
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2250,2500,3000]
stride=1
unit=0.005

stop_100 = 20000

nb_box=4

nbC = 32
nbO = 64

for V in Volumes

    for T in Temperatures

        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
        file="TRAJEC.xyz"

        if isfile( string(folder,file) )

            traj = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)
            nb_steps=size(traj)[1]
            nb_atoms=size(traj[1].names)[1]

            

        end
    end
end
