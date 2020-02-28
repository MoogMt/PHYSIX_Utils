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

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base="/media/mathieu/Elements/CO2/"


Volumes      = [ 10.0, 9.8, 9.5, 9.4, 9.375, 9.35, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05, 9.0, 8.82, 8.8, 8.6 ]
Temperatures = [ 1750, 2000, 2500, 3000 ]

for V in Volumes
    for T in Temperatures
        folder_target = string( folder_base, V, "/", T, "K/" )
        file_traj = string( folder_target, "/TRAJEC_fdb.xyz" )
        if ! isfile(file_traj)
            continue
        end
        print("V: ",V," T: ",T," K\n")
        folder_out=string(folder_target,"Data/Exp/")
        if ! isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end
        file_out = string(folder_out,"MSD.dat")
        msd = computingMSD( V, T, file_traj, file_out )
        file_outC = string( folder_out, "MSD_C.dat")
        msd_C =computingMSD( V, T, file_traj, file_outC, ["C"], [6] )
        file_outO = string( folder_out, "MSD_O.dat" )
        msd_O =computingMSD( V, T, file_traj, file_outO, ["O"], [8] )
        if size( msd_C )[1] == ()
            continue
        end
        msd_CO = msd_C-msd_O
        file_o = open( string(folder_out,"MSD_CO.dat"), "w" )
        for i=1:size( msd_C )[1]
            Base.write( file_o, string( i, " ", msd_CO[i], "\n" ) )
        end
        close( file_o )
    end
end
