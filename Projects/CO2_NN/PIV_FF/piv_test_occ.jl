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

nb_points=10000

test_matrix=zeros(nb_test,nb_train)
nb_piv_element=Int(nb_atoms*(nb_atoms-1)/2)

test_piv=zeros(nb_piv_element,nb_test)

d0=1.75
n=5

nbC=32
nbO=64


# Computing PIV_test
for test_step=1:nb_points
    print("Progress: ",test_step/nb_test*100,"%\n")
    # CC
    count_element = 1
    for carbon1=1:nbC-1
         for carbon2=carbon1+1:nbC
             @inbounds test_piv[count_element,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+start_point],cell,carbon1,carbon2) , d0, n )
            count_element += 1
        end
    end
    @inbounds test_piv[1:count_element,test_step]=sort(test_piv[1:count_element,test_step])
    # CO
    start_block=count_element-1
    count_element=1
    for carbon=1:nbC
        for oxygen=1:nbO
             @inbounds test_piv[count_element+start_block,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+start_point],cell,carbon,nbC+oxygen) , d0, n )
            count_element += 1
        end
    end
     @inbounds test_piv[start_block:start_block+count_element,test_step]=sort(test_piv[start_block:start_block+count_element,test_step])
    # OO
    start_block=count_element-1+start_block
    count_element=1
    for oxygen1=1:nbO-1
        for oxygen2=oxygen1+1:nbO
             @inbounds test_piv[count_element+start_block,test_step] = utils.switchingFunction( cell_mod.distance(traj[test_step+start_point],cell,nbC+oxygen1,nbC+oxygen2) , d0, n )
            count_element += 1
        end
    end
     @inbounds test_piv[start_block:start_block+count_element-1,test_step]=sort(test_piv[start_block:start_block+count_element-1,test_step])
end

# Computing test_matrix
for i=1:nb_point
    print("Progress: ",i/nb_test*100,"%\n")
    for j=1:nb_point
        dist=0
        for k=1:nb_piv_element
             @inbounds dist += (test_piv[k,i]-test_piv[k,j])*(test_piv[k,i]-test_piv[k,j])
        end
         @inbounds test_matrix[i,j]=sqrt(dist)
    end
end

file_test=open(string(folder_base,"matrix-",nb_point,"-d0_",d0,"-n_",n,"-m_",m,".dat"),"w")
for i=1:nb_test
    for j=1:nb_train
        Base.write(file_test,string(test_matrix[i,j]," "))
    end
    Base.write(file_test,"\n")
end
close(file_test)
