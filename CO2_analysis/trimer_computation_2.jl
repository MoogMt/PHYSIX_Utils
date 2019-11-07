GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

# Compute the number of occurences of trimers

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb

nbC=32
nbO=64

cut_off=1.75

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

Volumes=[9.325]
Temperatures=[1750,2000,2250,2500,2750,3000]
for T in Temperatures
    for V in Volumes

        folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
        folder_out=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/Data/")
        file=string(folder,"TRAJEC_wrapped.xyz")

        if isfile( string(folder,"TRAJEC_wrapped.xyz") )
            traj = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)
            nb_steps=size(traj)[1]
            nb_atoms=size(traj[1].names)[1]
            dimer_out=open(string(folder_out,"trimer-method2-",cut_off,".dat"),"w")
            for step=1:nb_steps
                print("Volume :",V," Temperature: ",T,"K Progress: ",step/nb_steps*100,"% \n")
                bonds=zeros(nbC,nbO)
                for carbon=1:nbC
                    for oxygen=1:nbO
                        if  cell_mod.distance( traj[step], cell, carbon, oxygen+nbC ) < cut_off
                            bonds[carbon,oxygen]=1
                        end
                    end
                end
                used=zeros(nbC)
                for carbon1=1:nbC-1
                    for carbon2=carbon1+1:nbC
                        for oxygen1=1:nbO
                            if bonds[carbon1,oxygen1] > 0  && bonds[carbon2,oxygen1] > 0
                                for carbon3=carbon2+1:nbC
                                    for oxygen2=1:nbO
                                        if bonds[carbon2,oxygen2] > 0  && bonds[carbon3,oxygen2] > 0
                                            for oxygen3=1:nbO
                                                if bonds[ carbon1, oxygen3] > 0 && bonds[carbon3, oxygen3] > 0
                                                    write(dimer_out,string(step," ",carbon1," ",carbon2," ",carbon3," ",oxygen1+nbC," ",oxygen2+nbC," ",oxygen3+nbC,"\n"))
                                                    break
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            close(dimer_out)
        end
    end
end
