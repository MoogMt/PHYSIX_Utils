using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod
using press_stress
using conversion
using correlation

using Statistics

function computeViscosity( stress_tensor::Array{T1,3}, cst::T2, max_lag::T3, dt::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Real  }
    corr_total=zeros(max_lag)
    count_=0
    for i=1:3
        for j=1:3
            if j != i

                for i=1:max_lag
                    corr_total[i] += corr[i]
                end
                count_ += 1
            end
        end
    end
    corr_total /= count_
    size_sig = size(corr)[1]
    int = 0
    for i=2:size_sig
        int += 0.5*(corr[i]+corr[i-1])*dt
    end
    return int*cst
end

kb = 1.38064852*10^(-23)
dt = 40*5*conversion.hatime2fs*10^(-12)

Volumes=[9.8,9.4,9.3,9.2,9.1,9.0,8.82]
Temperatures=[2000]

folder_base=string("/media/mathieu/Elements/CO2/")

folder_out = string( folder_base, "Data/Viscosity/" )
if ! isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end

# file_out = open( string(folder_out, "Visc_all.dat"), "w" )
# for T in Temperatures
#     file_t = open( string(folder_out, "Visc_",T,".dat"), "w" )
#     for V in Volumes

mc=12.0107/(5.489*10^(-4))
mo=15.999/(5.489*10^(-4))
V=8.82
T=3000
nbC=32
nbO=64

volume = V*V*V*conversion.a3tom3
folder_target = string( folder_base, V, "/", T, "K/" )
#cst = volume/(kb*T)
file_in = string( folder_target, "FTRAJECTORY_fdb" )
# if ! isfile( file_in )
#     continue
# end

x,v,f=cpmd.readFtraj( file_in )

file_out=open(string( folder_target, "/Data/Exp/Visc2.dat"),"w")
seg_length = 5000
nb_step = size(x)[1]
nb_seg=round(Int,nb_step/seg_length)
for step=1:seg_length
    visc_time = 0
    for seg=1:nb_seg
        seg_off = (seg-1)*seg_length
        value=0
        for carbon=1:nbC
            value += ( x[step+seg_off,carbon,1]*mc*v[step+seg_off,carbon,2]  ) - ( x[seg_off+1,carbon,1]*mc*v[seg_off+1,carbon,2]  )
        end
        for oxygen=1:nbO
            value += ( x[step+seg_off,nbC+oxygen,1]*mo*v[step+seg_off,nbC+oxygen,2]  ) - ( x[seg_off+1,nbC+oxygen,1]*mo*v[seg_off+1,nbC+oxygen,2]  )
        end
        value = value*value
        visc_time += value
    end
    write( file_out, string( step, " ", visc_time/nb_seg, "\n" ) )
end
close( file_out )
#         write( file_out, string( V, " ", T, " ", visc, "\n" ) )
#         write( file_t, string( V, " ", visc, "\n" ) )
#     end
#     close( file_t )
# end
# close( file_out )
