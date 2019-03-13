GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2500,3000]
Cut_Off=[1.75]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10
#

T=2000
V=9.8

# for V in Volumes
#     for T in Temperatures
# Folders
folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC.xyz")
folder_out=string(folder_in,"Data/")

    # if ! isfile(file)
    #     continue
    # end

# Reading xyz
print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]


nb_delta=Int(trunc(nb_steps/nb_cut))

MSD_C=zeros(nbC,nb_cut,nb_delta)
for carbon=1:nbC
    for cut=1:nb_cut
        start_cut=Int((cut-1)*nb_delta)+1
        end_cut=Int(cut*nb_delta)+1
        positions=traj[start_cut].positions[carbon,:]
        count=1
        for step=start_cut:end_cut-1
            dist=0
            for i=1:3
                dist += (traj[step].positions[carbon,i]-traj[start_cut].positions[carbon,i])*(traj[step].positions[carbon,i]-traj[start_cut].positions[carbon,i])
            end
            dist=sqrt(dist)
            MSD_C[carbon,cut,count] = dist*dist
            count+=1
        end
    end
end

MSD_C2=zeros(nb_delta)
for carbon=1:nbC
    for cut=1:nb_cut
        start_cut=Int((cut-1)*nb_delta)+1
        end_cut=Int(cut*nb_delta)+1
        positions=traj[start_cut].positions[carbon,:]
        count=1
        for step=start_cut:end_cut-1
            dist=0
            for i=1:3
                dist += (traj[step].positions[carbon,i]-traj[start_cut].positions[carbon,i])*(traj[step].positions[carbon,i]-traj[start_cut].positions[carbon,i])
            end
            dist=sqrt(dist)
            MSD_C[carbon,cut,count] = dist*dist
            count+=1
        end
    end
end


file_out=open(string(folder_out,"MSD_C_total.dat"),"w")
for step=1:size(MSD_C)[3]
    write(file_out,string(step," "))
    for carbon=1:size(MSD_C)[1]
        for cut=1:nb_cut
            write(file_out,string(MSD_C[carbon,cut,step]," "))
        end
    end
    write(file_out,string("\n"))
end
close(file_out)

file_out=open(string(folder_out,"MSD_C_Avg.dat"),"w")
for step=1:size(MSD_C)[3]
    write(file_out,string(step," "))
    MSD=0
    count=0
    for carbon=1:size(MSD_C)[1]
        for cut=1:nb_cut
            MSD += MSD_C[carbon,cut,step]
            count += 1
        end
    end
    write(file_out,string(" ",MSD/count,"\n"))
end
close(file_out)

MSD_O=zeros(nbO,nb_cut,nb_delta)
for oxygen=1:nbO
    for cut=1:nb_cut
        start_cut=Int((cut-1)*nb_delta)+1
        end_cut=Int(cut*nb_delta)+1
        positions=traj[start_cut].positions[nbC+oxygen,:]
        count=1
        for step=start_cut:end_cut-1
            dist=0
            for i=1:3
                dist += (traj[step].positions[nbC+oxygen,i]-traj[start_cut].positions[nbC+oxygen,i])*(traj[step].positions[nbC+oxygen,i]-traj[start_cut].positions[nbC+oxygen,i])
            end
            dist=sqrt(dist)
            MSD_O[oxygen,cut,count] = dist*dist
            count+=1
        end
    end
end

file_out=open(string(folder_out,"MSD_O_total.dat"),"w")
for step=1:size(MSD_O)[3]
    write(file_out,string(step," "))
    for oxygen=1:size(MSD_O)[1]
        for cut=1:nb_cut
            write(file_out,string(MSD_O[oxygen,cut,step]," "))
        end
    end
    write(file_out,string("\n"))
end
close(file_out)

file_out=open(string(folder_out,"MSD_O_Avg.dat"),"w")
for step=1:size(MSD_O)[3]
    write(file_out,string(step," "))
    MSD=0
    count=0
    for oxygen=1:size(MSD_O)[1]
        for cut=1:nb_cut
            MSD += MSD_O[oxygen,cut,step]
            count += 1
        end
    end
    write(file_out,string(" ",MSD/count,"\n"))
end
close(file_out)


file_out=open(string(folder_out,"MSD_CO_Avg.dat"),"w")
for step=1:size(MSD_C)[3]
    write(file_out,string(step," "))
    MSDC=0
    count=0
    for carbon=1:size(MSD_C)[1]
        for cut=1:nb_cut
            MSDC += MSD_C[carbon,cut,step]
            count += 1
        end
    end
    MSDC = MSDC/count
    count=0
    MSDO = 0
    for oxygen=1:size(MSD_O)[1]
        for cut=1:nb_cut
            MSDO += MSD_O[oxygen,cut,step]
            count += 1
        end
    end
    MSDO = MSDO/count
    write(file_out,string(" ",MSDO-MSDC,"\n"))
end
close(file_out)

# end
# end
