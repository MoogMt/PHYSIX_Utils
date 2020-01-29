# Loading file
include("contactmatrix.jl")

# Folder where to find file
folder="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/8.82/3000K/"

# number of centers
number_structure = 20

# Convergence parameters
n_iterations=1000000
k=0.02
s=0.005
kT=0.00000002

function computeEnergy( positions, real_distances, k )

    return energy
end

# Others
best_max_distance = 0.

# Reading frame to frame matrix
if isfile( string(folder,"FRAME_TO_FRAME.matrix") )

    # Reading file
    data=open(string(folder,"FRAME_TO_FRAME.matrix"))
    n, dmax = split(readline(data))
    n=Int(n)
    dmax=parse(Float64,dmax)
    d_real=zeros(n,n)
    for i=1:n
        print("Reading FRAME_TO_FRAME: ",i/n*100,"%\n")
        line=readline(data
        for j=1:n
            d_real[i,j]=parse(Float64,line[j])
        end
    end
    close(data)

    # Start with random positions
    pos=rand(n,2)

    # Compute and print Initial energy
    energy=0
    for i=1:n
        for j=1:n
            energy += 0.5*k*(sqrt((pos[i,1]-pos[j,1])^2+(pos[i,2]-pos[j,2])^2)-d_real[i,j])^2
        end
    end

    # Copy of position
    pos_bis=pos

    error_map=open(string(folder,"error_map.dat"),"w")
    conv=open(string(folder,"conv.dat"),"w")
    for iteration=1:n_iterations
        # Moving randomly a single element
        pos_bis[Int(trunc(rand()*n+1)),:] += (rand(2)-0.5)*s
        energy_bis=0
        max_distance=0
        for i=1:n
            for j=1:n
                d_map=sqrt((positions[i,1]-positions[j,1])^2+(positions[i,2]-positions[j,2])^2)
                energy_bis += 0.5*k*(d_map-d_real[i,j])^2
                if d_map > max_distance
                    max_distance = d_map
                end
                if iteration == n_iterations
                    print("error between ",i," and ",j," : ",(d_map-d_real[i,j]), ", relative: ",(d_map-d_real[i,j])/d_real[i,j],", orig_d mapped_d: ",d_map," ",d_real[i,j],"\n")
                    write(error_map,string(d_map," ",d_real[i,j],"\n")
                end
            end
        end
        # Evalutating change of energy
        de=(energy_bis-energy)/kT
        if de > 0
            if de < 50
                m=exp(-de)
            else
                m=0
            end
        else
            m=2.0
        end
        # Metropolis choice
        if rand() < m
            for i=1:n
                pos = pos_bis
                energy = energy_bis
            end
        end
        print("Progress: ",iteration/n_iterations*100,"%\n")
        write(conv,string(it," ",max_distance,"\n"))
    end
    close(error_map)
    close(conv)

    map=open(string(folder,",map.dat"),"w")
    for i=1:n
        for j=1:2
            write(map,string(pos[i,j]," "))
        end
        write(map,"\n")
    end
    close(map)
end
