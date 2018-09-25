include("contactmatrix.jl")

function sortmatrix( x )
    sizex=size(x)[1]
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                x[i]=x[j]
                x[j]=stock
            end
        end
    end
    return x
end

temperature=[2000,2250,2500,3000]
volume=[8.82,9.0,9.05,9.1,9.2,9.3,9.8]

cut_off=[1.6,1.7,1.8]

stride_sim=5
fs2ps=0.001
time_sim=1 # in fs
unit=time_sim*fs2ps*stride_sim# in ps
stride_analysis=1
start_time=0

for cutoff in cut_off
    for T in temperature
        for V in volume

            folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
            file_in=string(folder,"TRAJEC_wrapped.xyz")

            start_step=Int(start_time/(unit*stride_sim))

            atoms = filexyz.read( file_in, stride_analysis, start_step )
            cell  = cell_mod.Cell_param( V, V, V )

            nb_steps=size(atoms)[1]
            nb_atoms=size(atoms[1].names)[1]

            out=open(string(folder,"3CO_count_",V,"_",T,"_",cutoff,".dat"),"w")
            for step=1:nb_steps
                print("progress:",step/nb_steps*100,"%\n")
                count = 0
                for carbon=1:32
                    distances=zeros(nb_atoms-32)
                    for oxygen=33:nb_atoms
                        distances[oxygen-32]=cell_mod.distance(atoms[step],cell,carbon,oxygen)
                    end
                    if sortmatrix(distances)[3] < cutoff
                        count+=1
                    end
                end
                write(out,string(step*unit," ",count/32,"\n"))
            end
            close(out)

        end
    end
end
