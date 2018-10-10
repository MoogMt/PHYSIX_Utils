# Loading file
include("contactmatrix.jl")

nbC=32
nbO=64

cut_off=1.75

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

Volumes=[9.3]
Temperatures=[3000]

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
            trimer_out=open(string(folder_out,"trimer-",cut_off,".dat"),"w")
            for step=1:nb_steps
	    	print("V = ",V," T=",T,"K Progress=",step/nb_steps*100,"%\n")
                # Compute the adjacency matrix
                adjacency_matrix = zeros(nb_atoms,nb_atoms)
                for atom1=1:nb_atoms-1
                    for atom2=atom1+1:nb_atoms
                        if  cell_mod.distance( traj[step], cell, atom1, atom2 ) < cut_off
                            adjacency_matrix[ atom1, atom2 ] = 1
                            adjacency_matrix[ atom2, atom1 ] = 1
                        end
                    end
                end
                write(trimer_out,string(step," ",trace(adjacency_matrix*adjacency_matrix*adjacency_matrix)/6)," \n")
            end
            close(trimer_out)
        end
    end
end
