include("contactmatrix.jl")

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]
Cut_Off=[1.75]

# Time values
unit=0.005
frac=0.5

for T in Temperatures
    for V in Volumes

        folder_in=string(folder_base,V,"/",T,"K/")
        file=string(folder_in,"TRAJEC_wrapped.xyz")

        if ! isfile(folder_in,"TRAJEC_wrapped.xyz")
            continue
        end

        folder_out=string(folder_in,"Data/")

        atoms = filexyz.readFastFile(file)
        cell=cell_mod.Cell_param(V,V,V)

        nb_atoms=size(atoms[1].names)[1]
        nb_steps=size(atoms)[1]

        file_matrix=open(string(folder_out,"matrix.dat"),"w")
        for step=1:nb_steps
            print("V: ",V," T:",T,"K Progress:",step/nb_steps*100,"%\n")
            write(file_matrix,string(step*unit," "))
            for atom1=1:nb_atoms-1
                for atom2=atom1+1:nb_atoms
                    write(file_matrix,string(cell_mod.distance(atoms[step],cell,atom1,atom2)," "))
                end
            end
            write(file_matrix,string("\n"))
        end
        close(file_matrix)
    end
end
