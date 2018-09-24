# Loading file
include("contactmatrix.jl")

function sortIndex( indexes, x )
    sizex=size(x)[1]
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stock2=indexes[i]
                indexes[i]=indexes[j]
                x[i]=x[j]
                indexes[j]=stock2
                x[j]=stock
            end
        end
    end
    return
end


nbC=32
nbO=64

distances_C1=zeros(nbO)
distances_C2=zeros(nbO)
index1=zeros(Int,nbO)
index2=zeros(Int,nbO)

stop_100=20000

common_neighbor=zeros(Int,2)

Cut_Off=[1.6,1.7,1.75,1.8]
Volumes=[8.6,9.5]#[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[3000]#[2000,2250,2500,2750,3000]

for cut_off in Cut_Off
    for V in Volumes
        for T in Temperatures

            folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
            file=string(folder,"TRAJEC_wrapped.xyz")

            if isfile( string(folder,"TRAJEC_wrapped.xyz") )

                traj = filexyz.readFastFile(file)
                cell=cell_mod.Cell_param(V,V,V)

                nb_steps=size(traj)[1]
                nb_atoms=size(traj[1].names)[1]

                out_file=open(string(folder,"dimer_detect-",cut_off,".dat"),"w")

                for step=1:stop_100
                    if step > nb_steps
                        break
                    end
                    print("V: ",V," T: ",T," Progress: ",step/stop_100*100,"%\n")
                    count_dimer = 0
                    for carbon1 = 1:nbC

                        for oxygen = 1:nbO
                            distances_C1[oxygen] = cell_mod.distance(traj[step],cell,carbon1,oxygen)
                        end
                        for i=1:nbO
                            index1[i]=i
                        end
                        sortIndex(index1, distances_C1)
                        for carbon2 = carbon1+1:nbC
                            for oxygen = 1:nbO
                                distances_C2[oxygen] = cell_mod.distance(traj[step],cell,carbon2,oxygen)
                            end
                            for i=1:nbO
                                index2[i]=i
                            end
                            sortIndex(index2, distances_C2)
                            count = 0
                            for neighbor1=2:5
                                for neighbor2=2:5
                                    if index1[neighbor1] == index2[neighbor2]
                                        if cell_mod.distance(traj[step],cell,carbon1,index1[neighbor1]) < cut_off
                                            if  cell_mod.distance(traj[step],cell,carbon2,index2[neighbor2]) < cut_off
                                                count += 1
                                                common_neighbor[count] = index1[neighbor1]
                                            end
                                        end
                                    end
                                end
                            end
                            if count > 1
                                count_dimer += 1
                                write(out_file,string(step," "))
                                write(out_file,string(carbon1," ",carbon2," ",common_neighbor[1]," ",common_neighbor[2]," ")    )
                                write(out_file,string("\n"))
                            end
                        end
                    end
                    if count_dimer > 1
                        print("At step ",step," counted: ",count_dimer,"\n")
                    end
                end
                close(out_file)
            end
        end
    end
end
