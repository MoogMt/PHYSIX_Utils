include("contactmatrix.jl")

temperature=2000
volume=9.05
V=volume
T=temperature

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",volume,"/",temperature,"K/")
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

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

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

nn=zeros(nb_steps*32,5)
for step=1:nb_steps
    print("progress: ",step/nb_steps*100,"%\n")
    for carbon=1:32
        distances=zeros(nb_atoms-32)
        for oxygen=33:nb_atoms
            distances[oxygen-32] = cell_mod.distance(atoms[step],cell,carbon,oxygen)
        end
        distances=sortmatrix(distances)
        nn[step,1]=step
        for i=2:5
            nn[step*(carbon-1)+carbon,i]=distances[i-1]
        end
    end
end

min_dist=1.0
max_dist=3.0
nb_dist=150
ddist=(max_dist-min_dist)/nb_dist

folder="/home/moogmt/"

# Static
for k=1:4
    print("file ", k, " on 4\n")
    # Making histogram
    hist1D=zeros(nb_dist)
    for l=1:nb_steps*32
        for i=1:nb_dist
            if min_dist+(i-0.5)*ddist < nn[l,k+1] && min_dist+(i+0.5)*ddist > nn[l,k+1]
                hist1D[i] += 1
            end
        end
    end
    # Commit to file
    file=string("Static_NN-",k,"_V-",V,"_T-",T,".dat")
    check=open(string(folder,file),"w")
    for i=1:nb_dist
        write(check,string(i*ddist+min_dist," ",hist1D[i],"\n"))
    end
    close(check)
end

total_sim=nb_steps*unit
time_windows=[0.5,1,2,5,10]

for timewindow in time_windows

nb_window=Int(trunc(total_sim/timewindow)+1)

# Temporal
for k=1:4
    print("file ", k, " on 4\n")
    hist2d=zeros(nb_window,nb_dist)
    for l=1:nb_steps*32
        for i=1:nb_window
            if nn[l,1]*unit > (i-0.5)*timewindow && nn[l,1]*unit < (i+0.5)*timewindow
                for j=1:nb_dist
                    if min_dist+(j-0.5)*ddist < nn[l,k+1] && min_dist+(j+0.5)*ddist > nn[l,k+1]
                        hist2d[i,j] += 1
                    end
                end
            end
        end
    end
    file=string("Dyn_NN-",k,"_V-",V,"_T-",T,"_Window",timewindow,".dat")
    check=open(string(folder,file),"w")
    for i=1:nb_window
        for j=1:nb_dist
            write(check,string(i*timewindow," ",min_dist+j*ddist," ",hist2d[i,j],"\n"))
        end
        write(check,"\n")
    end
    close(check)
end

end
