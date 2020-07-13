# Loading file
using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using contact_matrix

nbC=32
nbO=64

folder_base = string( "/media/moogmt/Elements/CO2/" )

V=9.8
T=2000
run_nb=1

folder_target = string( folder_base, V, "/", T, "K/",run_nb,"-run/" )
file_traj = string( folder_target, "TRAJEC_db.xyz" )

traj = filexyz.readFileAtomList( file_traj )
cell = cell_mod.Cell_param( V, V, V )

nb_step = size(traj)[1]
nb_atoms = size(traj[1].names)[1]

for step=1:nb_step
    for atom=1:nb_atoms
        for i=1:3
            traj[step].positions[atom,i] = traj[step].positions[atom,i]%V
        end
    end
end

filexyz.writeXYZ( string(folder_target,"TRAJEC_test.xyz"),traj)

file_traj_wrap = string( folder_target, "TRAJEC_wrapped.xyz" )
traj_wrapped = filexyz.readFileAtomList( file_traj_wrap )

RMSD = zeros( nb_step )
for step=1:nb_step
    for atom=1:nb_atoms
        for i=1:3
            dist = traj_wrapped[step].positions[atom,i] - traj[step].positions[atom,i]
            RMSD[step] += dist*dist
        end
    end
    RMSD[step] = sqrt(RMSD[step])/nb_atoms
end

file_out = open(string("/home/moogmt/test.dat"),"w")
for step=1:nb_step
    write(file_out, string(step," ",RMSD[step],"\n"))
end
close(file_out)
