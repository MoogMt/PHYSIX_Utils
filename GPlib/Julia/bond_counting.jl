include("contactmatrix.jl")

temperature=3000
volume=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.375,9.4,9.8]

T=temperature

for V in volume

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file_in=string(folder,"TRAJEC_wrapped.xyz")

stride_sim=5
fs2ps=0.001
time_sim=0.5 # in fs
unit=time_sim*fs2ps*stride_sim# in ps

stride_analysis=1
start_time=5
start_step=Int(start_time/(unit*stride_sim))

atoms = filexyz.read( file_in, stride_analysis, start_step )
cell=cell_mod.Cell_param( V, V, V )

cut_off = [1.6,1.7,1.8]

for co in cut_off

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

nb_bonds=zeros(nb_steps)

out=open(string("/home/moogmt/bonds_count_",V,"_",T,"_",co,".dat"),"w")
for step=1:nb_steps
    print("progress : ",step/nb_steps*100,"%\n")
    count=0
    for carbon=1:32
        for oxygen=33:96
            if cell_mod.distance(atoms[step],cell,carbon,oxygen) < co
                count += 1
            end
        end
    end
    nb_bonds[step]=count
    write(out,string(step*unit," ",count,"\n"))
end
close(out)

total_sim=nb_steps*unit
time_windows=[0.5,1,2,5,10]

for timewindow in time_windows

    nb_window=Int(trunc(total_sim/timewindow)+1)

    hist1D=zeros(nb_window)

    for step=1:nb_steps
        for i=1:nb_window
            if step*unit > (i-0.5)*timewindow && step*unit < (i+0.5)*timewindow
                hist1D[i]+= nb_bonds[step]
            end
        end
    end

    out_box=open(string("/home/moogmt/bonds_count_",V,"_",T,"_",co,"A_box-",timewindow,".dat"),"w")
    for i=1:nb_window
        write(out_box,string(i*timewindow," ",hist1D[i]*unit/timewindow,"\n"))
    end
    close(out_box)

end
end
end
