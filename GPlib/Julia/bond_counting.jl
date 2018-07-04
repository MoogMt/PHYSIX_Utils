include("contactmatrix.jl")

temperature=3000
volume=9.05

V=volume
T=temperature

folder=string("/home/moogmt/CO2/CO2_AIMD/",volume,"/",temperature,"K/")
file_in=string(folder,"TRAJEC_wrapped.xyz")

stride_sim=5
fs2ps=0.001
time_sim=0.5 # in fs
unit=time_sim*fs2ps*stride_sim# in ps

stride_analysis=1
start_time=5
start_step=Int(start_time/(unit*stride_sim))

atoms = filexyz.read( file_in, stride_analysis, start_step )
cell=cell_mod.Cell_param( volume, volume, volume )

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
    write(out,string(step*unit," ",count,"\n"))
end
close(out)

end
