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

    traj = filexyz.readFileAtomList( file_traj )
    if traj == false
        return false
    end

    cell=cell_mod.Cell_param(V,V,V)

    msd_global=exp_data.computeMSD( traj, ["C","O"],[6.0,8.0] )

    return msd_global
end
function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString }

    msd  = computingMSD(V,T,file_traj)
    if msd == false
        return false
    end

    nb_step=size(msd)[1]
    file_o=open(file_out,"w")
    for step=1:nb_step
        Base.write(file_o,string(step," ",msd[step],"\n"))
    end
    close(file_o)

    return msd
end
function computingMSD( V::T1, T::T2, file_traj::T3, names::Vector{T4}, masses::Vector{T5} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: Real }

    traj = filexyz.readFileAtomList( file_traj )
    if traj == false
        return false
    end

    cell=cell_mod.Cell_param(V,V,V)

    msd_global=exp_data.computeMSD( traj, names, masses )

    return msd_global
end
function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4, names::Vector{T5}, masses::Vector{T6} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString, T6 <: Real }

    msd = computingMSD( V, T, file_traj, names, masses )
    if msd == false
        return false
    end

    nb_step=size(msd)[1]
    file_o=open(file_out,"w")
    for step=1:nb_step
        Base.write(file_o,string(step," ",msd[step],"\n"))
    end
    close(file_o)

    return msd
end

# Model for a linear fit
@. model(x, p) = p[1]*x

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base="/media/mathieu/Elements/CO2/"


Volumes=[9.8]
Temperatures=[3000]

for V in Volumes
    for T in Temperatures
        folder_target = string( folder_base, V, "/", T, "K/" )
        file_traj = string( folder_target, "/TRAJEC_fdb.xyz" )
        if ! isfile(file_traj)
            continue
        end
        folder_out=string(folder_target,"Data/Exp/")
        if ! isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end
        file_out = string(folder_out,"MSD.dat")
        msd,test = computingMSD( V, T, file_traj, file_out )
        file_outC = string( folder_out, "MSD_C.dat")
        msd_C,test=computingMSD( V, T, file_traj, file_outC, ["C"], [6] )
        file_outO = string( folder_out, "MSD_O.dat" )
        msd_O,test=computingMSD( V, T, file_traj, file_outO, ["O"], [8] )
        msd_CO = msd_C - msd_O
        file_o = open( string(folder_out,"MSD_CO.dat"), "w" )
        for i=1:size( msd_C )[1]
            Base.write( file_o, string( i, " ", msd_CO[i], "\n" ) )
        end
        close( file_o )
    end
end

# function blockAverage( data::Vector{T1}, max_block_size::T2 ) where { T1 <: Real, T2 <: Int }
#     size_record=[]
#     sigma_record=[]
#     for size_block=2:1:max_block_size
#         nb_block=Int(trunc(size(data)[1]/size_block))
#         meta_avg=0
#         avg_all_blocks=0
#         count_=1
#         for block=1:nb_block
#             avg_block=0
#             var_block=0
#             for i=1:size_block
#                 avg_block += data[count_]
#                 var_block += data[count_]*data[count_]
#                 count_ += 1
#             end
#             avg_block=avg_block/size_block
#             meta_avg+=sqrt(var_block/size_block - avg_block*avg_block)
#             avg_all_blocks = avg_all_blocks + avg_block
#         end
#         avg_all_blocks = avg_all_blocks/nb_block
#         sigma = sqrt( meta_avg/nb_block )
#         push!(size_record,size_block)
#         push!(sigma_record,sigma)
#     end
#     return size_record,sigma_record
# end
