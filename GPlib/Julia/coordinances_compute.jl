# Loading file
include("contactmatrix.jl")

unit=0.005
cut_off=1.75

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8,10.0]
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
            matrix=zeros(nb_atoms,nb_atoms)
            C2=0
            C3=0
            C4=0
            file_coord_local=open(string(folder_out,"Coordinances_time.dat"),"w")
            C2p=0
            C3p=0
            C4p=0
            for step=1:nb_steps
                matrix=zeros(nb_atoms,nb_atoms)
                print("V=",V," T=",T,"K Progress: ",step/nb_steps*100,"%\n")
                for i=1:nb_atoms-1
                    for j=i+1:nb_atoms
                        if cell_mod.distance(traj[step],cell,i,j) < cut_off
                            matrix[i,j] += 1
                            matrix[j,i] += 1
                        end
                    end
                end
                for i=1:nb_atoms
                    sum_i=sum(matrix[i,:])
                    if sum_i == 2
                        C2 += 1
                        C2p += 1
                    elseif sum_i == 3
                        C3 += 1
                        C3p += 1
                    elseif sum_i == 4
                        C4 += 1
                        C4p += 1
                    end
                end
                sum_local=C2p+C3p+C4p
                write(file_coord_local,string(step*unit," ",C2p/sum_local," ",C3p/sum_local," ",C4p/sum_local,"\n"))
            end
            close(file_coord_local)

            sum_all=C2+C3+C4
            file_coord=open(string(folder_out,"Avg_Coordinances.dat"),"w")
            write(file_coord,string(C2/sum_all," ",C3/sum_all," ",C4/sum_all,"\n"))
            close(file_coord)

        end
    end
end
