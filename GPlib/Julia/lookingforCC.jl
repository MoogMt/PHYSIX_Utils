# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=1
co=cut_off[1]

nbC=32
nbO=2*nbC

for V in volumes

    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
    file="TRAJEC_wrapped.xyz"
    print(folder,file,"\n")

    if isfile( string(folder,file) )
        traj = filexyz.read(string(folder,file),stride)
        cell=cell_mod.Cell_param(V,V,V)
        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]
        file=open(string("/home/moogmt/lookingCC-",V,"-",T,"-",co,".dat"),"w")
        for step=1:nb_steps
            for carbon1=1:nbC
                for carbon2=carbon1+1:nbC
                    if cell_mod.distance(traj[step],cell,carbon1,carbon2) < co
                        write(file,string(step," ",carbon1," ",carbon2," ",cell_mod.distance(traj[step],cell,carbon1,carbon2),"\n"))
                    end
                end
            end
        end
        close(file)
    else
        print("file: ",string(folder,file)," does not exists!\n")
    end

end
