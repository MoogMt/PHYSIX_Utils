# Loading file
include("contactmatrix.jl")

unit=0.005
cut_off=1.75

nbC=32
nbO=64

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
            O1=0
            O2=0
            file_coord_localC=open(string(folder_out,"CoordinancesC_time.dat"),"w")
            file_coord_localO=open(string(folder_out,"CoordinancesO_time.dat"),"w")

            for step=1:nb_steps
                C2p=0
                C3p=0
                C4p=0
                O1p=0
                O2p=0
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
                for i=1:nbC
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
                sum_localC=C2p+C3p+C4p
                for i=1:nbO
                    sum_i=sum(matrix[nbC+i,:])
                    if sum_i == 1
                        O1 += 1
                        O1p += 1
                    elseif sum_i == 2
                        O2 += 1
                        O2p += 1
                    end
                end
                sum_localO=O1p+O2p
                write(file_coord_localC,string(step*unit," ",C2p/sum_localC," ",C3p/sum_localC," ",C4p/sum_localO,"\n"))
                write(file_coord_localO,string(step*unit," ",O1p/sum_localO," ",O2p/sum_localO,"\n"))
            end
            close(file_coord_localC)
            close(file_coord_localO)

            sum_allC=C2+C3+C4
            file_coordC=open(string(folder_out,"Avg_CoordinancesC.dat"),"w")
            write(file_coordC,string(C2/sum_allC," ",C3/sum_allC," ",C4/sum_allC,"\n"))
            close(file_coordC)

            sum_allO=O1+O2
            file_coordO=open(string(folder_out,"Avg_CoordinancesO.dat"),"w")
            write(file_coordO,string(O1/sum_allO," ",O2/sum_allO,"\n"))
            close(file_coordO)
        end
    end
end
