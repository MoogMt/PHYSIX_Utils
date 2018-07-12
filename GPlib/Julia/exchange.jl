include("contactmatrix.jl")

temperature=3000
volume=[8.82]

T=temperature
#for V in volume
V=volume[1]
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

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

function sortmatrixwithindex( x )
    sizex=size(x)[1]
    indexes=zeros(sizex)
    for i=1:sizex
        indexes[i]=i
    end
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stock2=indexes[i]
                x[i]=x[j]
                indexes[i]=indexes[j]
                x[j]=stock
                indexes[j]=stock2
            end
        end
    end
    return x,indexes
end

nbC=32
nbO=64
cut_off=1.6

for oxygen=1:nbO
    print("Progress - Oxygen: ",oxygen/nbO*100,"%\n")
    neighbours=zeros(nb_steps,2)
    indexs=zeros(nb_steps,2)
    used=zeros(nb_steps)
    for step=1:nb_steps
        distances=zeros(32)
        for carbon=1:32
            distances[carbon] = cell_mod.distance(atoms[step],cell,carbon,oxygen+nbC)
        end
        distances, indexes = sortmatrixwithindex(distances)
        for i=1:2
            if distances[i] < cut_off
                neighbours[step,i] = 1
            else
                used[i] = 1
            end
        end
        indexs[step,:]=[indexes[1], indexes[2]]
    end
    count=0
    lifes=[]
    # Merge
    for neigh=1:2
        for i=step:nb_steps
            for j=1:nb_steps
            end
        end
    end
    # Stat
end

# Stat
