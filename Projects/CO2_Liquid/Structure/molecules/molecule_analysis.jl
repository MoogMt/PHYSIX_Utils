# Loading necessary stuff
using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix
using filexyz
using graph
using exp_data
using geom


# Thermodynamical values
Volumes=[10.0,9.8,9.5,9.4,9.375,9.35,9.325,9.3,9.25,9.2,9.15,9.1,9.05,9.0,8.82,8.8,8.6]
Temperature=[1750,2000,2500,3000,3500]
cut_off=1.75

V=8.82
T=3000


folder_base=string("/media/mathieu/Elements/CO2/",V,"/",T,"K/")

file_traj = string(folder_base,"TRAJEC_fdb_wrapped.xyz")

traj = filexyz.readFileAtomList( file_traj )
cell = cell_mod.Cell_param(V,V,V)
cut_off=1.75

nb_step=size(traj)[1]

for step=1:nb_step

    # Computing the molecules through graph exploration
    positions_local=copy(traj[step].positions)
    matrix = contact_matrix.buildMatrix( traj[step], cell, cut_off )
    molecules=graph.getGroupsFromMatrix(matrix)
    nb_molecules=size(molecules)[1]
    matrices=graph.extractAllMatrixForTrees( matrix, molecules )

    # For each molecule, compute add its size to the general histogram
    # And for the histogram of the infinite/finite size chains
    # Finally, prints unwrapped molecules in separate .xyz file
    # depending on if they are infinite or not.
    for molecule=1:nb_molecules
        if size(molecules[molecule])[1] <= 1
            continue
        end
        check = checkInfiniteChain( matrices[molecule], positions_local, cell, molecules[molecule], cut_off  )
        if check

        else

        end
    end

    traj[step].positions=positions_local
end
