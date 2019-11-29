GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]


nbC=32
nbO=nbC*2

V=9.8
T=3000

Volumes=[9.8]
Temperatures=[3000]

print("V=",V," T=",T,"\n")

# file_out_map=open(string(folder_base,"map_poisson.dat"),"w")
# for T in Temperatures
#     for V in Volumes
V=9.8
T=3000

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")
folder_out=string(folder_in,"Data/")

# if  ! isfile(file)
#     continue
# end

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

state_target=[0, 0, 0, 0, 0, 2, 1, 1, 0, 0]
target_number=-1
for i=1:size(states[1])[1]
    found=true
    for j=1:size(states[1])[2]
        if states[1][i,j] != state_target[j]
            found=false
            break
        end
    end
    if found
        global target_number=i
        break
    end
end

delta=200 # 200 steps = 1ps
d_delta=500
occurences_nb=[]
for step_start=delta:d_delta:nb_steps-delta
    print("Progress: ",step_start/(nb_steps-delta)*100,"%\n")
    for carbon=1:nbC
        occurence=0
        for step=step_start:step_start+delta
            if state_matrices[1][carbon,step] == target_number
                occurence += 1
                # Looking up next valid step
                for next=step+1:step_start+delta
                    if state_matrices[1][carbon,next] != target_number
                        check=true
                        for check_filter=next+1:next+10
                            if check_filter > nb_steps
                                break
                            end
                            if state_matrices[1][carbon,check_filter] == target_number
                                found2=false
                                for check_filter2=check_filter+1:step_start+delta
                                    if state_matrices[1][carbon,check_filter2] != target_number
                                        step=check_filter2
                                        found2=true
                                        break
                                    end
                                end
                                if found2
                                    check=false
                                    break
                                end
                            end
                        end
                        if check
                            step = next-1 # it will get incremeted at end of loop
                        end
                    end
                end
            end
        end
        push!(occurences_nb,occurence)
    end
end


file_out=open(string(folder_out,string("poisson-",delta,".dat")),"w")
for i=1:size(occurences_nb)[1]
    write(file_out,string(i," ",occurences_nb[i],"\n"))
end
close(file_out)

nb_=size(occurences_nb)[1]

lambda=sum(occurences_nb)/nb_

lambda2=0
count_=0
for i=1:nb_
    if occurences_nb[i] > 0
        global lambda2 += occurences_nb[i]
        global count_ += 1
    end
end
lambda2 /= count_

max_=0
for occ=1:nb_
    if max_ < occurences_nb[occ]
        global max_ = occurences_nb[occ]
    end
end

hist1D=zeros(Real,Int(max_)+1)
for occ=1:nb_
    hist1D[occurences_nb[occ]+1] += 1
end
hist1D[1]=exp(-lambda2)
hist1D/=sum(hist1D)

file_out=open(string(folder_out,string("poisson-hist-",delta,".dat")),"w")
for i=1:Int(max_)
    if i < 20
        write(file_out,string(i-1," ",hist1D[i]," ",exp(-lambda)*lambda^(i-1)/factorial(i-1)," ",exp(-lambda2)*lambda2^(i-1)/factorial(i-1),"\n"))
    else
        write(file_out,string(i-1," ",hist1D[i]," ",0," ",0,"\n"))
    end
end
close(file_out)

file_in_p=open(string(folder_out,"Avg_Pressure-BootStrap-nboot_1000.dat"))
lines=readlines(file_in_p)
close(file_in_p)

P=parse(Float64,split(lines[1])[2])

# write(file_out_map,string(P," ",T," ",lambda,"\n"))

#     end
# end
# close(file_out_map)


# unit=0.005
#
# lengths=[]
# for carbon=1:nbC
#     counting=false
#     count_=0
#     for step=1:nb_steps
#         if ! counting
#             # Found event
#             if state_matrices[1][carbon,step] != 1
#                 # We move at the end of the chain
#                 check=false
#                 for step_2=1:nb_steps
#                     if state_matrices[1][carbon,step_2] == 1 && step_2-step > 10
#                         step=step_2
#                         check=true
#                         break
#                     end
#                 end
#                 # If we did not find it, we go to the next carbon
#                 if ! check
#                     step=nb_steps+1
#                 end
#                 # And we start counting to the next event
#                 counting = true
#                 count_=1
#             end
#         else
#             # Counting steps to next event
#             if state_matrices[1][carbon,step] == 1
#                 count_ += 1
#             else
#                 if count_ > 1
#                     push!(lengths,count_)
#                 end
#                 counting=false
#                 count_=0
#             end
#         end
#     end
# end

# min_value=unit
# max_value=unit
# for i=1:size(lengths)[1]
#     if max_value < lengths[i]
#         global max_value = lengths[i]
#     end
# end
#
# # Histogram
# nb_box=50
# delta=(max_value-min_value)/nb_box
# hist1D=zeros(Real,nb_box)
# for i=1:size(lengths)[1]
#     hist1D[ Int(trunc( lengths[i]/delta-min_value ))+1 ] += 1
# end
# hist1D/=sum(hist1D)
#
# # Writting data
# file_out=open(string(folder_out,"hist_poisson.dat"),"w")
# for i=1:nb_box
#     write(file_out,string( (i*delta+min_value)*unit," ",hist1D[i],"\n"))
# end
# close(file_out)




#end
