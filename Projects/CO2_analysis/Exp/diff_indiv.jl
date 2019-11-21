GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using exp_data
using LsqFit

function computingMSD( V::T1, T::T2, file_traj::T3 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString }

    traj,test=filexyz.readFastFile(file_traj)
    if ! test
        return zeros(1,1), false
    end

    cell=cell_mod.Cell_param(V,V,V)

    msd_global=exp_data.computeMSD( traj, ["C","O"],[6.0,8.0] )

    return msd_global, true
end
function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString }

    msd,test=computingMSD(V,T,file_traj)
    if ! test
        return zeros(1,1), false
    end

    nb_step=size(msd)[1]
    file_o=open(file_out,"w")
    for step=1:nb_step
        Base.write(file_o,string(step," ",msd[step],"\n"))
    end
    close(file_o)

    return msd, true
end
function computingMSD( V::T1, T::T2, file_traj::T3, names::Vector{T4}, masses::Vector{T5} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: Real }

    traj,test=filexyz.readFastFile(file_traj)
    if ! test
        return zeros(1,1), false
    end

    cell=cell_mod.Cell_param(V,V,V)

    msd_global=exp_data.computeMSD( traj, names, masses )

    return msd_global, true
end
function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4, names::Vector{T5}, masses::Vector{T6} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString, T6 <: Real }

    msd,test=computingMSD( V, T, file_traj, names, masses )
    if ! test
        return zeros(1,1), false
    end

    nb_step=size(msd)[1]
    file_o=open(file_out,"w")
    for step=1:nb_step
        Base.write(file_o,string(step," ",msd[step],"\n"))
    end
    close(file_o)

    return msd, true
end

# Model for a linear fit
@. model(x, p) = p[1]*x

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"

Volumes=[8.82]
Temperatures=[3000]

for V in Volumes
    for T in Temperatures
        file_traj=string(folder_base,"/",V,"/",T,"K/1-run/TRAJEC.xyz")
        file_out=string(folder_base,"/",V,"/",T,"K/Data/MSD.dat")
        msd,test=computingMSD(V,T,file_traj,file_out)
        file_outC=string(folder_base,"/",V,"/",T,"K/Data/MSD_C.dat")
        msd_C,test=computingMSD(V,T,file_traj,file_outC,["C"],[6])
        file_outO=string(folder_base,"/",V,"/",T,"K/Data/MSD_O.dat")
        msd_O,test=computingMSD(V,T,file_traj,file_outO,["O"],[8])
    end
end


# D_avg=[]
# for start_cut=1:nb_space:nb_steps-nb_delta
#     print("Progress",start_cut/(nb_steps-nb_delta)*100,"%\n")
#     MSD_local=zeros(nb_delta)
#     for atom=1:nb_atoms
#         positions=traj[start_cut].positions[atom,:]
#         for step=1:nb_delta
#             dist=0
#             for i=1:3
#                 dist += ( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )*( (traj[step+start_cut].positions[atom,i]-traj[start_cut].positions[atom,i]) - (barycenter_all[step+start_cut,i]-barycenter_all[start_cut,i]) )
#             end
#             MSD_local[step] += dist
#         end
#     end
#     MSD_local /= (nbC+nbO)
#
#     times_MSD=clustering.simpleSequence(nb_delta)*0.005
#
#     fit = curve_fit(model, times_MSD , MSD_local , [0.4] )
#     push!(D_avg,coef(fit)[1])
# end


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
