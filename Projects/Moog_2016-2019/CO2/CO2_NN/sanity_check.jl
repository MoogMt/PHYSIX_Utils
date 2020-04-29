using utils
using atom_mod
using cell_mod
using filexyz
using pdb
using conversion
using cpmd
using press_stress
using exp_data


folder_base=string("/media/moogmt/Elements/CO2/")

V=9.8
T=2000
folder_target = string( folder_base, V, "/", T, "K/1-run/")

file_traj = string( folder_target, "TRAJEC_db.xyz" )
file_ftraj = string( folder_target, "FTRAJECTORY_db")

positions,velocities,forces=cpmd.readFtraj( file_ftraj )
traj = filexyz.readFileAtomList( file_traj )

positions *= conversion.bohr2Ang

cell = cell_mod.Cell_param(V,V,V)

nb_atoms=96

target_step=1
for atom=1:nb_atoms
    for i=1:3
        positions[target_step,atom,i] = cell_mod.wrap( positions[target_step,atom,i],V)
    end
end
