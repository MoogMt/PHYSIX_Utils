GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"clustering.jl"))


using LsqFit

@. model(x, p) = p[1]*x

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]
Temperatures=[2000,2500,3000]

# Cut-off distance for bonds
cut_off_bond = 1.75
nbC=32
nbO=64
nb_cut=10

V=8.82
T=3000

# for T in Temperatures
#     for V in Volumes
folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC.xyz")
folder_out=string(folder_in,"Data/")

print(V," ",T,"\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_steps=size(traj)[1]
nb_atoms=size(traj[1].names)[1]

barycenter_C=zeros(nb_steps,3)
barycenter_O=zeros(nb_steps,3)
barycenter_all=zeros(nb_steps,3)
for step=1:nb_steps
    for carbon=1:nbC
        for i=1:3
            barycenter_C[step,i] += traj[step].positions[carbon,i]
            barycenter_all[step,i] += traj[step].positions[carbon,i]
        end
    end
    for oxygen=1:nbO
        for i=1:3
            barycenter_O[step,i] += traj[step].positions[nbC+oxygen,i]
            barycenter_all[step,i] += traj[step].positions[nbC+oxygen,i]
        end
    end
    for j=1:3
        barycenter_C[step,j] = barycenter_C[step,j]/nbC
        barycenter_O[step,j] = barycenter_O[step,j]/nbO
        barycenter_all[step,j] = barycenter_all[step,j]/(nbC+nbO)
    end
end


nb_delta=5000
nb_space=50

MSD_local=zeros(nb_delta)

D_avg=[]
for start_cut=1:nb_space:nb_steps-nb_delta
    print("Progress",start_cut/(nb_steps-nb_delta)*100,"%\n")
    MSD_local=zeros(nb_delta)
    for atom=1:nb_atoms
        positions=traj[start_cut].positions[atom,:]
        for step=1:nb_delta
            dist=0
            for i=1:3
                dist += ( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )*( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )
            end
            MSD_local[step] += dist
        end
    end
    MSD_local /= (nbC+nbO)

    times_MSD=clustering.simpleSequence(nb_delta)*0.005

    fit = curve_fit(model, times_MSD , MSD_local , [0.4] )
    push!(D_avg,coef(fit)[1])
end

D2=zeros(size(D_avg)[1])
for i=1:size(D_avg)[1]
    D2[i]=D_avg[i]
end

function blockAverage( data::Vector{T1}, max_block_size::T2 ) where { T1 <: Real, T2 <: Int }
    size_record=[]
    sigma_record=[]
    for size_block=2:1:max_block_size
        nb_block=Int(trunc(size(data)[1]/size_block))
        meta_avg=0
        avg_all_blocks=0
        count_=1
        for block=1:nb_block
            avg_block=0
            var_block=0
            for i=1:size_block
                avg_block += data[count_]
                var_block += data[count_]*data[count_]
                count_ += 1
            end
            avg_block=avg_block/size_block
            meta_avg+=sqrt(var_block/size_block - avg_block*avg_block)
            avg_all_blocks = avg_all_blocks + avg_block
        end
        avg_all_blocks = avg_all_blocks/nb_block
        sigma = sqrt( meta_avg/nb_block )
        push!(size_record,size_block)
        push!(sigma_record,sigma)
    end
    return size_record,sigma_record
end

size_,sigma_=blockAverage(D2,140)

file_out=open(string(folder_out,"D_",nb_delta,"-",nb_space,"-block.dat"),"w")
for i=1:size(sigma_)[1]
    write(file_out,string(size_[i]," ",sigma_[i],"\n"))
end
close(file_out)

# end
# end
