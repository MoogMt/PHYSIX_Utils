path_to_libatomsim = "C:\\Users\\moogm\\OneDrive\\Documents\\GitHub\\LibAtomicSim\\Julia\\"
push!( LOAD_PATH, path_to_libatomsim )

using atom_mod
using cell_mod
using geom
using cube_mod
using periodicTable
using conversion

using LinearAlgebra

path_to_trimer_folder = "F:\\PHYSIX\\Mathieu\\CO2\\AIMD\\Liquid\\PBE-MT\\ELF\\ELF\\Rings\\"

frame_num=50

path_to_file = string( path_to_trimer_folder, "ELF_", frame_num, "_treated.cube")

atoms, cell, volume = cube_mod.readCube( path_to_file )

output_path_file = string( path_to_trimer_folder, "ELF_", frame_num, "_treated_mod.cube" )

cube_mod.writeCube( output_path_file, atoms, cell, volume )
