using filexyz
using cpmd
using utils
using contact_matrix
using atom_mod
using cell_mod
using press_stress
using conversion

using Statistics

# Compute the pressure

Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperatures=[3000,2500,2000,1750]

folder_base=string("/media/mathieu/Elements/CO2/")

nb_atoms=96
masses=[6,8]
nb_species=[32,64]

total_mass=0
for i=1:2
    global total_mass += masses[i]*nb_species[i]
end
total_mass *= conversion.amu2gram

for T in Temperatures
    folder_out_general = string( folder_base, "Data/Pressure/" )
    if !isdir( folder_out_general )
        Base.Filesystem.mkdir(folder_out_general )
    end
    file_out_general = string( folder_out_general,T,"_P_global.dat")
    handle_out_gen = open( file_out_general, "w" )
    for V in Volumes
        folder_target=string(folder_base, V, "/",T,"K/")
        file_target = string(folder_target,"Pressure_fdb")
        if ! isfile(file_target)
            continue
        end
        pressure = press_stress.readPressure( file_target )
        folder_out = string( folder_target, "Data/Pressure/")
        if !isdir( folder_out )
            Base.Filesystem.mkdir(folder_out)
        end
        file_out_loc = string( folder_out, "P_global.dat" )
        handle_out_loc = open( file_out_loc, "w")
        Base.write( handle_out_loc, string( Statistics.mean(pressure)," ",Statistics.std(pressure),"\n"))
        Base.write( handle_out_gen, string( V," ", V*V*V/(nb_atoms), " ", Statistics.mean(pressure)," ",Statistics.std(pressure)," ",total_mass/(V*V*V*conversion.a3tocm3),"\n"))
        close( handle_out_loc )
    end
    close( handle_out_gen )
end
