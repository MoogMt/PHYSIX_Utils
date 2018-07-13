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
oxygens2=[]
carbons2=[]
starts=[]
end_at_endsim=[]

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
            if neighbours[i,1] > 0
                neighbours[i,1] = 0
                print("Progress - Oxygen: ",oxygen/nbO*100," ")
                print("Counting Life:",i/nb_steps*100,"%\n")
                j=i+1
                push!(list_index,indexs[i,neigh])
                push!(list_start,i)
                while j <= nb_steps
                    if indexs[j,1] == indexs[i,neigh] && neighbours[j,1] > 0
                        neighbours[j,1] = 0
                        j+=1
                    elseif indexs[j,2] == indexs[i,neigh] && neighbours[j,1] > 0
                        neighbours[j,2] = 0
                        j+=1
                    else
                        break
                    end
                end
                push!(list_end,j)
            end
        end
    end
    used=zeros(size(list_index)[1])
    for i=1:size(list_index)[1]
        if used[i] == 0
            used[i] = 1
            for j=i+1:size(list_index)[1]
                if used[j] == 0
                    if list_index[i] == list_index[j]
                        list_end[i]=list_end[j]
                        used[j]=1
                    end
                end
            end
            push!(lifes,list_end[i]-list_start[i])
            push!(oxygens2,oxygen)
            push!(carbons2,list_index[i])
            push!(starts,list_start[i])
        end
    end
end

file=open(string(),"w")
for i=1:size(oxygens2)[1]
    write(file,string("O:",oxygens2[i]," C:",carbons2[i]," life:",lifes[i]," or ",lifes[i]*unit,"ps, start at step:",starts[i]+start_step))
    if end_at_endsim[i] == 0
        write(file,string())
    else
        write(file,string())
    end
    write(file,"\n")
end
close(file)
