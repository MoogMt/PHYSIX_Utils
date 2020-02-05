GPfolder=string("/home/mathieu/LibAtomicSim/Julia/")
push!(LOAD_PATH, GPfolder)
B2O3folder=string("/home/mathieu/B2O3_Transitions/")
push!(LOAD_PATH,B2O3folder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using cpmd
using resfile
using pdb
using filexyz
using utils
using topomap
using pimaim

using b2o3database

# Basic informations
#------------------------------------------------------------------------------
folder_base=string("/home/mathieu/B2O3/")
structure_folder_path = string(folder_base,"Structures/GuillaumeFerlat/")
#------------------------------------------------------------------------------

files = utils.getAllFilesWithExtension( structure_folder_path, "pdb" )
names=utils.getFilesName(files)

file_path = string( structure_folder_path , "B2O3-I/B2O3-I_720.pdb")

if ! isfile(file_path)
    print("OUPS ! ")
else
    atoms,cell = pdb.readStructureAtomList( file_path )
    cell = cell_mod.nonOrtho2Ortho( cell_mod.params2Matrix( cell ) )
    pdb.writePdb( string(structure_folder_path,"f22/f22_init.pdb"), atoms, cell )
    pimaim.writeRestart( string(structure_folder_path,"f22/f22_init_restart.dat"), atoms, cell_mod.params2Matrix(cell) )
    pimaim.writeCrystalCell( string(structure_folder_path,"f22/f22_init_crystal_cell.inpt"), atoms, cell_mod.params2Matrix(cell) )
end
