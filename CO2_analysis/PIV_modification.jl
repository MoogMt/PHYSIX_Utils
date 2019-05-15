GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"utils.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/CO2/CO2_AIMD/"


T=3000
V=8.82

n_run=1

folder_in=string(folder_base,V,"/",T,"K/",n_run,"-run/")
folder_in=string(folder_base,V,"/",T,"K/")

print("Computing Data\n")
traj=filexyz.readFastFile(string(folder_in,"TRAJEC.xyz"))
cell=cell_mod.Cell_param(V,V,V)

traj=traj[1:5]
nb_structure=size(traj)[1]
nb_atoms=size(traj[1].names)[1]
nbC=32
nbO=64

# Wrapping
for i=1:nb_structure
    cell_mod.wrap(traj[i],cell)
end

d0=1.75
n=5

# ComputePIV
piv_element=Int(nb_atoms*(nb_atoms-1)/2)
piv=zeros(piv_element,nb_structure)
for step=1:nb_structure
    print("PIV computation progress: ",step/nb_structure*100,"%\n")
    count=1
    start=count
    for carbon=1:nbC
        for carbon2=carbon+1:nbC
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,carbon2),d0,n)
            count = count + 1
        end
    end
    piv[start:count-1,step]=sort(piv[start:count-1,step])
    start=count
    for oxygen=1:nbO
        for oxygen2=oxygen+1:nbO
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,nbC+oxygen,nbC+oxygen2),d0,n)
            count = count + 1
        end
    end
    # Sort
    piv[start:count-1,step]=sort(piv[start:count-1,step])
    start=count
    for carbon=1:nbC
        for oxygen=1:nbO
            piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,nbC+oxygen),d0,n)
            count = count + 1
        end
    end
    # Sort
    piv[start:count-1,step]=sort(piv[start:count-1,step])
end

# Naive PIV
dist1=0
for i=1:piv_element
    global dist1 +=  ( piv[i,1]-piv[i,2])*(piv[i,1]-piv[i,2])
end

nb_box=200
delta=1/nb_box
piv2=zeros(piv_element,nb_structure)
for step=1:2
    hist=zeros(Int,nb_box)
    for carbon=1:nbC-1
        for carbon2=carbon+1:nbC
            hist[Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,carbon2),d0,n)/delta))+1] += 1
        end
    end
    # Reconstruct PIV
    count_local=0
    for i=1:nb_box
        piv2[1+count_local:count_local+hist[i],step] = ones(Int,hist[i])*i
        count_local += hist[i]
    end
    hist=zeros(Int,nb_box)
    for oxygen=nbC+1:nbC+nbO-1
        for oxygen2=oxygen+1:nbC+nbO
            hist[Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,oxygen,oxygen2),d0,n)/delta))+1] += 1
        end
    end
    count_box=Int(nbC*(nbC-1)/2)
    count_local=0
    for i=1:nb_box
        piv2[count_local+count_box+1:count_local+count_box+hist[i],step] = ones(hist[i])*i
        count_local += hist[i]
    end
    hist=zeros(Int,nb_box)
    for carbon=1:nbC-1
        for oxygen=nbC+1:nbC+nbO
            hist[Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,oxygen),d0,n)/delta))+1] += 1
        end
    end
    count_local=0
    count_box=Int(nbC*(nbC-1)/2+nbO*(nbO-1)/2)
    for i=1:nb_box
        piv2[count_local+count_box+1:count_local+count_box+hist[i],step] = ones(hist[i])*i
        count_local += hist[i]
    end
end

# Actual PIV
dist2=0
div_=1/(delta*delta)
for i=1:3
    global dist2 +=  ( piv[i,1]-piv[i,2])*(piv[i,1]-piv[i,2])*div_
end


nb_box=200
hist=zeros(2,nb_box,3)
delta=1/nb_box
for step=1:2
    for carbon=1:nbC-1
        for carbon2=carbon+1:nbC
            hist[step,Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,carbon2),d0,n)/delta))+1,1] += 1
        end
    end
    for oxygen=nbC+1:nbC+nbO-1
        for oxygen2=oxygen+1:nbC+nbO
            hist[step,Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,oxygen,oxygen2),d0,n)/delta))+1,2] += 1
        end
    end
    for carbon=1:nbC-1
        for oxygen=nbC+1:nbC+nbO
            hist[step,Int(trunc(utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,oxygen),d0,n)/delta))+1,3] += 1
        end
    end
end

# Distance 1-2
dist3=0
div_=(delta*delta)
hist_cal=zeros(nb_box)
for k=1:3
    delta=hist[1,:,k]-hist[2,:,k]
    for i=1:nb_box
        if delta[i] > 0
            for j=i+1:nb_box
                if delta[i] == 0
                    break
                end
                if hist[2,j,k] == 0
                    continue
                end
                if hist[2,j,k] >= abs(delta[i])
                    hist[2,j,k] = hist[2,j,k] - delta[i]
                    delta[j] = delta[j] + delta[i]
                    global dist3 = dist3 + delta[i]*delta[i]*(j-i)*div_
                    print("add: ", delta[i]*delta[i]*(j-i)*div_,"\n")
                    break
                else
                    print("add2: ",hist[2,j,k]*hist[2,j,k]*(j-i)*div_," ",hist[2,j,k],"\n")
                    global dist3 = dist3 + hist[2,j,k]*hist[2,j,k]*(j-i)*div_
                    hist[2,j,k] = 0
                end
            end
        elseif delta[i] < 0
            for j=i+1:nb_box
                if delta[i] == 0
                    break
                end
                if hist[1,j,k] == 0
                    continue
                end
                if hist[1,j,k] >= abs(delta[i])
                    hist[1,j,k] = hist[1,j,k] - delta[i]
                    delta[j] = delta[j] + delta[i]
                    global dist3 = dist3 + delta[i]*delta[i]*(j-i)*div_
                    print("add3: ", delta[i]*delta[i]*(j-i)*div_,"\n")
                    break
                else
                    global dist3 = dist3 + hist[1,j,k]*hist[1,j,k]*(j-i)*div_
                    print("add2: ",hist[1,j,k]*hist[1,j,k]*(j-i)*div_,"\n")
                    hist[1,j,k] = 0
                end
                if delta[i] == 0
                    break
                end
            end
        end
    end
end



# Naive way to compute forces
# distances=zeros(nb_structure,nb_structure)
# for step1=1:nb_structure
#     print("Computing distance matrix: ",step1/nb_structure*100,"%\n")
#     for step2=step1+1:nb_structure
#         for i=1:size(piv)[1]
#             distances[step1,step2] += (piv[i,step1]-piv[i,step2])*(piv[i,step1]-piv[i,step2])
#         end
#         distances[step1,step2] = sqrt(distances[step1,step2])
#         distances[step2,step1] = distances[step1,step2]
#     end
# end
#
#
# for step=1:nb_structure
#     print("PIV computation progress: ",step/nb_structure*100,"%\n")
#     count=1
#     start=count
#     for carbon=1:nbC
#         for carbon2=carbon+1:nbC
#             piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,carbon2),d0,n)
#             count = count + 1
#         end
#     end
#     piv[start:count-1,step]=sort(piv[start:count-1,step])
#     start=count
#     for oxygen=1:nbO
#         for oxygen2=oxygen+1:nbO
#             piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,nbC+oxygen,nbC+oxygen2),d0,n)
#             count = count + 1
#         end
#     end
#     # Sort
#     piv[start:count-1,step]=sort(piv[start:count-1,step])
#     start=count
#     for carbon=1:nbC
#         for oxygen=1:nbO
#             piv[count,step] = utils.switchingFunction(cell_mod.distance(traj[step],cell,carbon,nbC+oxygen),d0,n)
#             count = count + 1
#         end
#     end
#     # Sort
#     piv[start:count-1,step]=sort(piv[start:count-1,step])
# end



#
# file_energy=open(string(folder_in,"EKS_base"))
# lines=readlines(file_energy)
# close(file_energy)
#
# energy=zeros(nb_structure)
# for i=1:nb_structure
#     energy[i] = parse(Float64,split(lines[i])[3])
# end
#
# file_out=open(string(folder_in,"PIV_1/results",d0,"-",n,".dat"),"w")
# for i=1:nb_structure
#     print("Progress: ",i/nb_structure*100,"%\n")
#     for j=i+1:nb_structure
#         write(file_out,string(distances[i,j]," ",abs(energy[i]-energy[j])/32*13.605693,"\n"))
#     end
# end
# close(file_out)
#
# file_out=open(string(folder_in,"PIV_1/piv-",d0,"-",n,".dat"),"w")
# for i=1:piv_element
#     print("Progress: ",i/piv_element*100,"%\n")
#     write(file_out,string(i," "))
#     for j=1:nb_structure
#         write(file_out,string(piv[i,j]," "))
#     end
#     write(file_out,string("\n"))
# end
# close(file_out)
