# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=5
co=cut_off[1]

nbC=32
nbO=2*nbC


    V=9.375
    T=3000
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
    file="TRAJEC_wrapped_nostride.xyz"

    if isfile( string(folder,file) )
        traj = filexyz.read(string(folder,file),stride)
        cell=cell_mod.Cell_param(V,V,V)
        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]
        file=open(string(folder,"traj_stride.xyz"),"w")
        for step=1:nb_steps
            write(file,string(nb_atoms,"\n"))
            write(file,string("CHECK\n"))
            for i=1:nb_atoms
                write(file,string(traj[step].names[i]," ",))
                for j=1:3
                    write(file,string(traj[step].positions[i,j]," ",))
                end
                write(file,string("\n"))
            end
        end
        close(file)
    else
        print("file: ",string(folder,file)," does not exists!\n")
    end
