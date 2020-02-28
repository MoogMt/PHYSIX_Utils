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
    corr = correlation.autocorr( stress_tensor[ :, 1, 2 ], max_lag )
    size_sig = size(corr)[1]
    int = 0
    for i=2:size_sig
        int += 0.5*(corr[i]+corr[i-1])*dt
    end
    return int*cst
end

kb = 1.38064852*10^(-23)
max_lag=10000
dt = 40*5*conversion.hatime2fs*10^(-12)

Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperatures=[3000,2500,2000,1750]

folder_base=string("/media/mathieu/Elements/CO2/")

folder_out = string( folder_base, "Data/Viscosity/" )
if ! isdir( folder_out )
    Base.Filesystem.mkdir( folder_out )
end

file_out = open( string(folder_out, "Visc_all.dat"), "w" )
for T in Temperatures
    file_t = open( string(folder_out, "Visc_",T,".dat"), "w" )
    for V in Volumes
        volume = V*V*V*conversion.a3tom3
        folder_target = string( folder_base, V, "/", T, "K/" )
        cst = volume/(kb*T)
        stress_tensor = readStress( string( folder_target, "STRESS_fdb" ) )
        visc = computeViscosity( stress_tensor, cst, max_lag, dt )
        write( file_out, string( V, " ", T, " ", visc, "\n" ) )
        write( file_t, string( V, " ", visc, "\n" ) )
    end
    close( file_t )
end
close( file_out )
