# Loading file
include("contactmatrix.jl")

func="PBE-MT"

Temperatures=[ 2000, 2250, 25000, 3000 ]
Volumes = [ 8.6, 8.82, 9.0, 9.05, 9.1, 9.15, 9.2, 9.25, 9.3, 9.325, 9.35, 9.375, 9.4, 9.5, 9.8 ]
Cut_Off = [ 1.6 , 1.7 , 1.75, 1.8 ]

nbC = 32
nbO = 64

unit=0.005

for cut_off in Cut_Off
    for V in Volumes
        for T in Temperatures

            folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"/")
            file=string(folder,"TRAJEC_wrapped.xyz")

            if ! isfile(string(folder,file))
                break
            end

            atoms = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            nb_steps=size(atoms)[1]
            nb_atoms=size(atoms[1].names)[1]

            file_C = open(string(folder,"C_coordinance.dat"),"w")
            file_O = open(string(folder,"O_coordinance.dat"),"w")

            for step=1:nb_steps

                print("Progress: ",step/nb_steps*100," %\n")

                write(file_C,string(step*unit," "))
                coordC_avg=0
                for carbon=1:nbC
                    coord_C=0
                    for atom=1:nb_atoms
                        if carbon != atom
                            if dist=cell_mod.distance(atoms[k],cell,carbon,atom) < cut_off
                                coord_C += 1
                            end
                        end
                    end
                    write(file_C,string(coord_C," "))
                    coordC_avg += coord_C
                end
                write(file_C,string(coordC_avg/nbC,"\n"))


                write(file_O,string(step*unit," "))
                coordO_avg=0
                for oxygen=1:nbO
                    coord_O = 0
                    for atom=1:nb_atoms
                        if oxygen+nbC != atom
                            if  dist=cell_mod.distance(atoms[k],cell,carbon,atom) < cut_off
                                coord_O += 1
                            end
                        end
                    end
                    write(file_O,string(coord_O," "))
                    coordC_avg += coord_O
                end
                write(file_C,string(coordO_avg/nbO,"\n"))

            end

            close(file_C)
            close(file_O)

        end
    end
end
