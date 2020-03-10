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
using Statistics
using exp_data


function readPressure( file_pressure::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_pressure )
        return false, false, false
    end
    handle_press_in = open( file_pressure )
    nb_p = 0
    while !eof( handle_press_in )
        readline( handle_press_in )
        nb_p += 1
    end
    seekstart( handle_press_in )
    pressure = zeros(Real, nb_p )
    delta_p  = zeros(Real, nb_p )
    volumes  = zeros(Real, nb_p )
    for i=1:nb_p
        keyword=split(readline( handle_press_in ))
        pressure[i] = parse( Float64, keyword[3] )
        delta_p[i]  = parse( Float64, keyword[4] )
        volumes[i]  = parse( Float64, keyword[1] )
    end
    close( handle_press_in )
    return pressure, delta_p, volumes
end

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"
folder_base="/media/moogmt/Elements/CO2/"
folder_base="/media/mathieu/Elements/CO2/"

nbC=32
nbO=64
size_sim = 20000 # in nb_step

Volumes      = [ 9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.1,9.05, 9.0, 8.82, 8.6 ]
Temperatures = [ 2000, 2500, 3000 ]

#10.0, 9.8, 9.5, 9.4, 9.375, 935, 9.325, 9.3, 9.25, 9.2, 9.15, 9.1, 9.05,
@. model(x, p) = p[1]+p[2]*x
p0 = [0.5, 0.5]

timestep = 0.005 # ps

d_conv = 10^(-8)*10^(-8)/(100*10^(-12))

block_size=2000  # 2000 steps = 10ps
n_block = round(Int, size_sim/block_size )
times_ = zeros(block_size)
for step=1:block_size
    times_[step] = (step-1)*timestep
end
folder_out_map = string( folder_base, "Data/Exp/" )
if ! isdir( folder_out_map )
    Base.Filesystem.mkdir( folder_out_map )
end
for T in Temperatures
    file_pressure = string( folder_base, "Data/Pressure/",T,"_P_global.dat")
    P,dP,Volumes = readPressure( file_pressure )
    handle_C_map = open( string( folder_out_map, "D_C_map_", T, "-", block_size, ".dat" ), "w" )
    handle_O_map = open( string( folder_out_map, "D_O_map_", T, "-", block_size, ".dat" ), "w" )
    handle_all_map = open( string( folder_out_map, "D_all_map_", T, "-", block_size, ".dat" ), "w" )
    for i_vol=1:size(Volumes)[1]
        folder_target = string( folder_base, Volumes[i_vol], "/", T, "K/" )
        file_traj = string( folder_target, "/TRAJEC_fdb.xyz" )
        if ! isfile( file_traj )
            continue
        end
        print( "V: ",Volumes[i_vol]," T: ",T,"K"," block_size: ",block_size,"\n")
        folder_out=string(folder_target,"Data/Exp/")
        if ! isdir( folder_out )
            Base.Filesystem.mkdir( folder_out )
        end
        traj = filexyz.readFileAtomList( file_traj )
        msd_total = zeros( n_block, block_size )
        msd_C = zeros( n_block, block_size )
        msd_O = zeros( n_block, block_size)
        for i_block=1:n_block
            traj_loc = traj[ (i_block-1)*block_size+1:i_block*block_size ]
            positions = atom_mod.getPositionsAsArray( traj_loc )
            bar_all = exp_data.computeBarycenter( positions, traj[1].names, ["C","O"], [6,8] )
            bar_C = exp_data.computeBarycenter( positions )
            bar_O = exp_data.computeBarycenter( positions )
            msd_total[ i_block, : ] = exp_data.computeMSD( positions, bar_all )
            msd_C[ i_block, : ] = exp_data.computeMSD( positions[ :, 1:nbC, : ], bar_C )
            msd_O[ i_block, : ] = exp_data.computeMSD( positions[ :, nbC+1:nbC+nbO, : ], bar_O )
        end
        file_out_all = string( folder_out, "MSD_all_", block_size, ".dat" )
        file_out_C   = string( folder_out, "MSD_C_",   block_size, ".dat" )
        file_out_O   = string( folder_out, "MSD_O_",   block_size, ".dat" )
        file_out_CO   = string( folder_out, "MSD_CO_",   block_size, ".dat" )
        handle_all_out = open( file_out_all, "w" )
        handle_C_out = open( file_out_C, "w" )
        handle_O_out = open( file_out_O, "w" )
        handle_CO_out = open( file_out_CO, "w")
        # msd_all2 = zeros( block_size )
        # msd_C2 = zeros( block_size )
        # msd_O2 = zeros( block_size )
        for step=1:block_size
            Base.write( handle_all_out, string( step*timestep, " ", Statistics.mean( msd_total[ :, step ] ), " ", Statistics.std( msd_total[ :, step ] )/sqrt(n_block), "\n" ) )
            Base.write( handle_C_out, string( step*timestep, " ", Statistics.mean( msd_C[ :, step ] ), " ", Statistics.std( msd_C[ :, step ] )/sqrt( n_block ), "\n" ) )
            Base.write( handle_O_out, string( step*timestep, " ", Statistics.mean( msd_O[ :, step ] ), " ", Statistics.std( msd_O[ :, step ] )/sqrt( n_block ), "\n" ) )
            Base.write( handle_CO_out, string( step*timestep, " ", Statistics.mean( msd_C[ :, step ] ) - Statistics.mean( msd_O[ :, step ] ) , "\n" ) )
            # msd_all2[ step ] = Statistics.mean( msd_total[ :, step ] )
            # msd_C2[ step ] = Statistics.mean( msd_C[ :, step ] )
            # msd_O2[ step ] = Statistics.mean( msd_O[ :, step ] )
        end
        close( handle_all_out )
        close( handle_C_out )
        close( handle_O_out )
        close( handle_CO_out )
        # fit = curve_fit(model, times_, msd_all2, p0 )
        # d_all = fit.param[2]
        # fit = curve_fit(model, times_, msd_C2, p0 )
        # d_coefC = fit.param[2]
        # fit = curve_fit(model, times_, msd_O2, p0 )
        # d_coefO = fit.param[2]
        d_all=zeros( n_block )
        d_coefC=zeros( n_block )
        d_coefO=zeros( n_block )
        for i_bloc=1:n_block
            fit = curve_fit(model, times_, msd_total[i_bloc,:], p0 )
            d_all[i_bloc] = fit.param[2]
            fit = curve_fit(model, times_, msd_C[i_bloc,:], p0 )
            d_coefC[i_bloc] = fit.param[2]
            fit = curve_fit(model, times_, msd_O[i_bloc,:], p0 )
            d_coefO[i_bloc] = fit.param[2]
        end
        # Base.write( handle_C_map, string( P[i_vol], " ", dP[i_vol], " ", d_all, "\n" ) )
        # Base.write( handle_O_map, string( P[i_vol], " ", dP[i_vol], " ", d_coefC, "\n" ) )
        # Base.write( handle_all_map, string( P[i_vol], " ", dP[i_vol], " ", d_coefO, "\n" ) )
        Base.write( handle_C_map, string( P[i_vol], " ", dP[i_vol], " ", Statistics.mean(d_all)*d_conv, " ", Statistics.std(d_all)*d_conv, "\n" ) )
        Base.write( handle_O_map, string( P[i_vol], " ", dP[i_vol], " ", Statistics.mean(d_coefC)*d_conv, " ", Statistics.std(d_coefC)*d_conv, "\n" ) )
        Base.write( handle_all_map, string( P[i_vol], " ", dP[i_vol], " ", Statistics.mean(d_coefO)*d_conv, " ", Statistics.std(d_coefO)*d_conv,  "\n" ) )
    end
    close( handle_C_map )
    close( handle_O_map )
    close( handle_all_map )
end


# function computingMSD( V::T1, T::T2, file_traj::T3 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString }
#
#     traj = filexyz.readFileAtomList( file_traj )
#     if traj == false
#         return false
#     end
#
#     cell=cell_mod.Cell_param(V,V,V)
#
#     msd_global=exp_data.computeMSD( traj, ["C","O"],[6.0,8.0] )
#
#     return msd_global
# end
# function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString }
#
#     msd  = computingMSD(V,T,file_traj)
#     if msd == false
#         return false
#     end
#
#     nb_step=size(msd)[1]
#     file_o=open(file_out,"w")
#     for step=1:nb_step
#         Base.write(file_o,string(step," ",msd[step],"\n"))
#     end
#     close(file_o)
#
#     return msd
# end
# function computingMSD( V::T1, T::T2, file_traj::T3, names::Vector{T4}, masses::Vector{T5} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: Real }
#
#     traj = filexyz.readFileAtomList( file_traj )
#     if traj == false
#         return false
#     end
#
#     cell=cell_mod.Cell_param(V,V,V)
#
#     msd_global=exp_data.computeMSD( traj, names, masses )
#
#     return msd_global
# end
# function computingMSD( V::T1, T::T2, file_traj::T3, file_out::T4, names::Vector{T5}, masses::Vector{T6} ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString, T6 <: Real }
#
#     msd = computingMSD( V, T, file_traj, names, masses )
#     if msd == false
#         return false
#     end
#
#     nb_step=size(msd)[1]
#     file_o=open(file_out,"w")
#     for step=1:nb_step
#         Base.write(file_o,string(step," ",msd[step],"\n"))
#     end
#     close(file_o)
#
#     return msd
# end

# for V in Volumes
#     for T in Temperatures
#         folder_target = string( folder_base, V, "/", T, "K/" )
#         file_traj = string( folder_target, "/TRAJEC_fdb.xyz" )
#         if ! isfile(file_traj)
#             continue
#         end
#         print("V: ",V," T: ",T," K\n")
#         folder_out=string(folder_target,"Data/Exp/")
#         if ! isdir( folder_out )
#             Base.Filesystem.mkdir( folder_out )
#         end
#         file_out = string(folder_out,"MSD.dat")
#         msd = computingMSD( V, T, file_traj, file_out )
#         file_outC = string( folder_out, "MSD_C.dat")
#         msd_C =computingMSD( V, T, file_traj, file_outC, ["C"], [6] )
#         file_outO = string( folder_out, "MSD_O.dat" )
#         msd_O =computingMSD( V, T, file_traj, file_outO, ["O"], [8] )
#         if size( msd_C )[1] == ()
#             continue
#         end
#         msd_CO = msd_C-msd_O
#         file_o = open( string(folder_out,"MSD_CO.dat"), "w" )
#         for i=1:size( msd_C )[1]
#             Base.write( file_o, string( i, " ", msd_CO[i], "\n" ) )
#         end
#         close( file_o )
#     end
# end
