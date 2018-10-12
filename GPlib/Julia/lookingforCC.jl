# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[9.325]
temperatures=[2000,2500,3000]
cut_off=[1.75]

func=functionnals[1]
stride=1

nbC=32
nbO=2*nbC

stop_step=20000

for co in cut_off
    for T in temperatures
        for V in volumes

            folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
            folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/Data/")
            file="TRAJEC_wrapped.xyz"
            print(folder,file,"\n")

            if isfile( string(folder,file) )
                traj = filexyz.read(string(folder,file),stride)
                cell=cell_mod.Cell_param(V,V,V)
                nb_steps=size(traj)[1]
                nb_atoms=size(traj[1].names)[1]
                file=open(string(folder_out,"CCbond-",V,"-",T,"-",co,".dat"),"w")
                for step=1:stop_step
                    if step > nb_steps
                        break
                    end
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
    end
end
