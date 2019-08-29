GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using utils
using geom
using atom_mod
using cell_mod
using filexyz
using cube_mod
using pdb

# Folder for data
V=8.82
T=3000

folder_base=string("/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/1-run/")
file="TRAJEC_wrapped.xyz"

# Getting Positions
traj = filexyz.readFastFile( string(folder_base,file))
cell = cell_mod.Cell_param( V, V, V )

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

start_point=2000

nb_train=8000
nb_test=1000

test_matrix=zeros(nb_test,nb_train)

nb_piv_element=Int(nb_atoms*(nb_atoms-1)/2)

train_piv=zeros(nb_piv_element,nb_train)
test_piv=zeros(nb_piv_element,nb_test)

d0=1.75
n=5

nbC=32
nbO=64

# Computing PIV_train
for train_step=1:nb_train
    print("Progress: ",train_step/nb_train*100,"%\n")
    # CC
    count_element = 1
    for carbon1=1:nbC-1
        for carbon2=carbon1+1:nbC
            @inbounds train_piv[count_element,train_step] = utils.switchingFunction( cell_mod.distance(traj[train_step+start_point],cell,carbon1,carbon2) , d0, n )
            count_element += 1
        end
    end
    train_piv[1:count_element,train_step]=sort(train_piv[1:count_element,train_step])
    # CO
    start_block=count_element-1
    count_element=1
    for carbon=1:nbC
        for oxygen=1:nbO
            @inbounds train_piv[count_element+start_block,train_step] = utils.switchingFunction( cell_mod.distance(traj[train_step+start_point],cell,carbon,nbC+oxygen) , d0, n )
            count_element += 1
        end
    end
    train_piv[start_block:start_block+count_element,train_step]=sort(train_piv[start_block:start_block+count_element,train_step])
    # OO
    start_block=count_element-1+start_block
    count_element=1
    for oxygen1=1:nbO-1
        for oxygen2=oxygen1+1:nbO
            @inbounds train_piv[count_element+start_block,train_step] = utils.switchingFunction( cell_mod.distance(traj[train_step+start_point],cell,nbC+oxygen1,nbC+oxygen2) , d0, n )
            count_element += 1
        end
    end
    train_piv[start_block:start_block+count_element-1,train_step]=sort(train_piv[start_block:start_block+count_element-1,train_step])
end

# file_test=open(string(folder_base,"test_piv.dat"),"w")
# for i=1:nb_piv_element
#     Base.write(file_test,string(i," ",train_piv[i,1],"\n"))
# end
# close(file_test)

# Computing PIV_test
 for test_step=1:nb_test
    print("Progress: ",test_step/nb_test*100,"%\n")
    # CC
    count_element = 1
    for carbon1=1:nbC-1
         for carbon2=carbon1+1:nbC
             @inbounds test_piv[count_element,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+nb_train+1+start_point],cell,carbon1,carbon2) , d0, n )
            count_element += 1
        end
    end
    @inbounds test_piv[1:count_element,test_step]=sort(test_piv[1:count_element,test_step])
    # CO
    start_block=count_element-1
    count_element=1
    for carbon=1:nbC
        for oxygen=1:nbO
             @inbounds test_piv[count_element+start_block,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+nb_train+1+start_point],cell,carbon,nbC+oxygen) , d0, n )
            count_element += 1
        end
    end
     @inbounds test_piv[start_block:start_block+count_element,test_step]=sort(test_piv[start_block:start_block+count_element,test_step])
    # OO
    start_block=count_element-1+start_block
    count_element=1
    for oxygen1=1:nbO-1
        for oxygen2=oxygen1+1:nbO
             @inbounds test_piv[count_element+start_block,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+nb_train+1+start_point],cell,nbC+oxygen1,nbC+oxygen2) , d0, n )
            count_element += 1
        end
    end
     @inbounds test_piv[start_block:start_block+count_element-1,test_step]=sort(test_piv[start_block:start_block+count_element-1,test_step])
end

# Computing train matrix
# for i=1:nb_train
#     print("Progress: ",i/nb_train*100,"%\n")
#     Threads.@threads for j=i+1:nb_train
#         dist=0
#         for k=1:nb_piv_element
#             @inbounds dist+= (train_piv[k,i]-train_piv[k,j])*(train_piv[k,i]-train_piv[k,j])
#         end
#         @inbounds train_matrix[i,j]=sqrt(dist)
#         @inbounds train_matrix[j,i]=train_matrix[i,j]
#     end
# end
#
# file_train=open(string(folder_base,"train_matrix-",nb_train,".dat"),"w")
# for i=1:nb_train
#     for j=1:nb_train
#         Base.write(file_train,string(train_matrix[i,j]," "))
#     end
#     Base.write(file_train,"\n")
# end
# close(file_train)


# Computing test_matrix
for i=1:nb_test
    print("Progress: ",i/nb_test*100,"%\n")
    for j=1:nb_train
        dist=0
        for k=1:nb_piv_element
             @inbounds dist += (test_piv[k,i]-train_piv[k,j])*(test_piv[k,i]-train_piv[k,j])
        end
         @inbounds test_matrix[i,j]=sqrt(dist)
    end
end

file_test=open(string(folder_base,"test_matrix-",nb_test,".dat"),"w")
for i=1:nb_test
    for j=1:nb_test
        Base.write(file_test,string(test_matrix[i,j]," "))
    end
    Base.write(file_test,"\n")
end
close(file_test)


file_in=open(string(folder_base,"ENERGIES"))

stride_=5

for i=1:start_point
    for j=1:stride_
        readline(file_in)
    end
end

# Reading energies train
energy_train=zeros(nb_train)
for i=1:nb_train
    energy_train[i]=parse(Float64,split(readline(file_in))[4])
    for j=1:stride_-1
        readline(file_in)
    end
end

# Reading energies test
energy_test=zeros(nb_test)
for i=1:nb_test
    energy_test[i]=parse(Float64,split(readline(file_in))[4])
    for j=1:stride_-1
        readline(file_in)
    end
end
close(file_in)

energy_prediction=zeros(nb_test)

function sigfunction( distance::T1, param::T2 ) where {T1 <: Real, T2 <: Int }
    a=distance^param
    return 1/a
end

param=3
for i=1:nb_test
    sum_sig=0
    for j=1:nb_train
        sig=sigfunction(test_matrix[i,j],param)
        sum_sig += sig
        energy_prediction[i] += sig*energy_train[j]
    end
    energy_prediction[i] = energy_prediction[i]/sum_sig
end

using conversion

avg_test=0
var_test=0
min_test=0
for i=1:nb_test
    if min_test > energy_test[i]
        global min_test=energy_test[i]
    end
    global avg_test += energy_test[i]
    global var_test += energy_test[i]*energy_test[i]
end
avg_test /= nb_test
var_test = sqrt(var_test/nb_test - avg_test*avg_test)

# Checking results
err=0
file_out=open(string(folder_base,"test_ff-",nb_train,"-",nb_test,"-",param,".dat"),"w")
for i=1:nb_test
    Base.write(file_out,string(i," ",(energy_prediction[i]-avg_test)*13.605693009/nbC," ",(energy_test[i]-avg_test)*13.605693009/nbC,"\n"))
    global err += (energy_prediction[i]*13.605693009/nbC-energy_test[i]*13.605693009/nbC)*(energy_prediction[i]*13.605693009/nbC-energy_test[i]*13.605693009/nbC)
end
close(file_out)
print("erreur: ",err/nb_test*1000,"\n")
