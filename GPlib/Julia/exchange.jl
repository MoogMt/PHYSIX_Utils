include("contactmatrix.jl")

temperature=3000
volume=[9.8]

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

lifes=[]
count=0

for oxygen=1:nbO
    print("Progress - Oxygen: ",oxygen/nbO*100,"%\n")
    neighbours=zeros(nb_steps,2)
    indexs=zeros(nb_steps,2)
    for step=1:nb_steps
        distances=zeros(32)
        for carbon=1:32
            distances[carbon] = cell_mod.distance(atoms[step],cell,carbon,oxygen+nbC)
        end
        distances, indexes = sortmatrixwithindex(distances)
        for i=1:2
            if distances[i] < cut_off
                neighbours[step,i] = 1
            end
        end
        indexs[step,:]=[indexes[1], indexes[2]]
    end
    count=0
    list_index=[]
    list_start=[]
    list_end=[]
    # Aggregate
    for neigh=1:2
        for i=1:nb_steps
            print("Progress - Oxygen: ",oxygen/nbO*100," ")
            print("Counting Life:",i/nb_steps*100,"%\n")
            index_check=indexs[i,neigh]
            j=i+1
            push!(list_index,index_check)
            push!(list_start,i)
            while ( indexs[i,1] == index_check || indexs[i,2] == index_check ) && j <= nb_steps
                j+=1
            end
            push!(list_end,j)
        end
    end
    # Wrap
    used = zeros( size(list_index)[1] )
    for i=1:size(list_index)[1]
        if used[i] == 0
            used[i] = 1
            for j=i+1:size(list_index)[1]
                if list_index[i] == list_index[j] && used[j] == 0
                    if list_end[i] - list_start[j] < 10
                        used[j]=1
                        list_end[i]=list_end[j]
                    end
                end
            end
            count+=1
            push!(lifes, list_end[i]-list_start[i] )
        end
    end
    # Stat
end

# Stat
