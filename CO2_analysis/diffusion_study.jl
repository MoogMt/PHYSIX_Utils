GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"clustering.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

func="PBE-MT"

unit=0.005
nbC = 32
nbO = 64

V=9.8
T=2000

folder=string(folder_base,V,"/",T,"K/")
folder_out=string(folder,"Data/")
file="TRAJEC.xyz"

traj = filexyz.readFastFile( string(folder,file) )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size( traj )[1]
nb_atoms=size( traj[1].names )[1]


barycenter=zeros(Real,nb_steps,3)
C_barycenter=zeros(Real,nb_steps,3)
O_barycenter=zeros(Real,nb_steps,3)
d_barycenter=zeros(Real,nb_steps)
d_barycenterC=zeros(Real,nb_steps)
d_barycenterO=zeros(Real,nb_steps)

file_out=open(string(folder_out,"barycenter_diff.dat"),"w")
for step=1:nb_steps
    print("Progress:",step/nb_steps*100,"\n")
    for i=1:3
        for carbon=1:nbC
            barycenter[step,i] += traj[step].positions[carbon,i]
            C_barycenter[step,i] += traj[step].positions[carbon,i]
        end
        C_barycenter[step,i] /= nbC
        for oxygen=1:nbO
            barycenter[step,i] += traj[step].positions[nbC+oxygen,i]
            O_barycenter[step,i] += traj[step].positions[nbC+oxygen,i]
        end
        O_barycenter[step,i] /= nbO
        barycenter[step,i] /= nb_atoms
    end
    for i=1:3
        d_barycenter[step] += (barycenter[step,i]-barycenter[1,i])*(barycenter[step,i]-barycenter[1,i])
        d_barycenterC[step] += (C_barycenter[step,i]-C_barycenter[1,i])*(C_barycenter[step,i]-C_barycenter[1,i])
        d_barycenterO[step] += (O_barycenter[step,i]-O_barycenter[1,i])*(O_barycenter[step,i]-O_barycenter[1,i])
    end
    d_barycenter[step] = sqrt(d_barycenter[step])
    d_barycenterC[step] = sqrt(d_barycenterC[step])
    d_barycenterO[step] = sqrt(d_barycenterO[step])
    write(file_out,string(step," ",d_barycenter[step]," ",d_barycenterC[step]," ",d_barycenterO[step],"\n"))
end
close(file_out)

barycenter=[]
C_barycenter=[]
O_barycenter=[]

file_position_C=open(string(folder_out,"positions_C.dat"),"w")
for step=1:nb_steps
    write(file_position_C,string(step," "))
    for carbon=1:nbC
        for i=1:3
            write(file_position_C,string(traj[step].positions[carbon,i]," "))
        end
    end
    write(file_position_C,string("\n"))
end
close(file_position_C)

file_position_O=open(string(folder_out,"positions_O.dat"),"w")
for step=1:nb_steps
    write(file_position_O,string(step," "))
    for oxygen=nbC+1:nb_atoms
        for i=1:3
            write(file_position_O,string(traj[step].positions[oxygen,i]," "))
        end
    end
    write(file_position_O,string("\n"))
end
close(file_position_O)

cut_off_distance=0.4

nb_points=5000

rho=zeros(Real,nb_points)
for i=1:nb_points
    print("Progress: ",i/nb_points*100,"%\n")
    for j=i+1:nb_points
        distance=0
        for k=1:3
            distance += (traj[i].positions[1,k]-traj[j].positions[1,k])*(traj[i].positions[1,k]-traj[j].positions[1,k])
        end
        gauss = exp( - distance/(cut_off_distance*cut_off_distance))
        rho[i] += gauss
        rho[j] += gauss
    end
end

index=clustering.simpleSequence(size(rho)[1])
for i=1:size(rho)[1]
	for j=i+1:size(rho)[1]
		if rho[i] < rho[j]
			stock=rho[i]
			rho[i]=rho[j]
			rho[j]=stock
			stock=index[i]
			index[i]=index[j]
			index[j]=stock
		end
	end
end

max_distance=0
for i=1:nb_points
	print("Progress: ",i/nb_points*100,"%\n")
	for j=i+1:nb_points
		distance=0
		for k=1:3
			distance += (traj[i].positions[1,k]-traj[j].positions[1,k])*(traj[i].positions[1,k]-traj[j].positions[1,k])
		end
		if max_distance < distance
			global max_distance = distance
		end
	end
end

delta=ones(nb_points)*max_distance
delta[ 1 ] = -1
nneigh=zeros(Int,nb_points)
for i=2:nb_points
	print("Progress: ",i/nb_points*100,"%\n")
	for j=1:i-1
		distance=0
		for k=1:3
			distance += (traj[ index[i] ].positions[1,k]-traj[ index[j] ].positions[1,k])*(traj[ index[i] ].positions[1,k]-traj[ index[j] ].positions[1,k])
		end
		if distance < delta[ i ]
			delta[ i ] = distance
			nneigh[ i ] = j
		end
	end
end

delta[1] = 0
for i=1:size(delta)[1]
	if delta[i] > delta[1]
		delta[1] = delta[i]
	end
end

file_out=open(string(folder_out,"decision-",1,"-",cut_off_distance,"-",nb_points,".dat"),"w")
for i=1:size(rho)[1]
	write(file_out,string(rho[i]," ",delta[i],"\n"))
end
close(file_out)

min_delta=5
min_rho=5

n_cluster=0
cl=ones(Int,nb_points)*(-1)
icl=[]
# Determine the cluster centers
for i=1:nb_points
    if rho[i] > min_rho && delta[i] > min_delta
		global n_cluster += 1
        cl[index[i]] = n_cluster
        global icl=push!(icl,index[i])
    end
end

for i=1:nb_points
    if cl[index[i]] == -1
        cl[index[i]] = cl[ index[nneigh[i]]  ]
    end
end

for i=1:n_cluster
    file=open(string(folder_out,"coord_predict-",1,"-",cut_off_distance,"-",nb_points,"-cluster-",i,"-index.dat"),"w")
    for j=1:nb_points
        if i == cl[j]
            for k=1:3
                write(file,string( traj[j].positions[1,k]," ",))
            end
            write(file,string("\n"))
        end
    end
    close(file)
end
