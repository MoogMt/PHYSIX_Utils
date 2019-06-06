GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.6]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]


nbC=32
nbO=nbC*2

V=9.8
T=3000

print("V=",V," T=",T,"\n")


for T in [2000,2500,3000]

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

print("Reading Trajectory\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]


cut_off_bond = 1.75
max_neigh=5

data,types,type_atoms=buildCoordinationMatrix( traj , cell , cut_off_bond, max_neigh )
nb_types=size(types)[1]
states, state_matrices, counts = assignDataToStates( data, nb_types, type_atoms )

writeStates(string(folder_out,"statesC-",cut_off_bond,".dat"),states[1],counts[1],types)
writeStates(string(folder_out,"statesO-",cut_off_bond,".dat"),states[2],counts[2],types)
writeStateMatrix(string(folder_out,"statesC-matrix-",cut_off_bond,".dat"),state_matrices[1])
writeStateMatrix(string(folder_out,"statesO-matrix-",cut_off_bond,".dat"),state_matrices[2])

unit=0.005

lengths=[]
for carbon=1:nbC
    counting=false
    count_=0
    for step=1:nb_steps
        if ! counting
            # Found event
            if state_matrices[1][carbon,step] != 1
                # We move at the end of the chain
                check=false
                for step_2=1:nb_steps
                    if state_matrices[1][carbon,step_2] == 1
                        step=step_2
                        check=true
                        break
                    end
                end
                # If we did not find it, we go to the next carbon
                if ! check
                    step=nb_steps+1
                end
                # Anti-flickering measure:
                # We start counting ONLY if the event lasted more than 4 frames
                # so at least 20fs

                # And we start counting to the next event
                counting = true
                count_=1
            end
        else
            # Counting steps to next event
            if state_matrices[1][carbon,step] == 1
                count_ += 1
            else
                if count_ > 1
                    push!(lengths,count_)
                end
                counting=false
                count_=0
            end
        end
    end
end

min_value=unit
max_value=unit
for i=1:size(lengths)[1]
    if max_value < lengths[i]
        max_value = lengths[i]
    end
end

# Histogram
nb_box=200
delta=(max_value-min_value)/nb_box
hist1D=zeros(Real,nb_box)
for i=1:size(lengths)[1]
    hist1D[ Int(trunc( lengths[i]/nb_box ))+1 ] += 1
end
hist1D/=sum(hist1D)

# Writting data
file_out=open(string(folder_out,"hist_poisson.dat"),"w")
for i=1:nb_box
    write(file_out,string( (i*delta+min_value)*unit," ",hist1D[i],"\n"))
end
close(file_out)

end
