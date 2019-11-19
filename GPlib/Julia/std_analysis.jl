module std_analysis

using atom_mod
using cell_mod

function computeRMSD( structure1::T1, structure2::T2 ) where { T1 <: atom_mod.AtomList, T2 <: AtomList }
    nb_atoms1=size(structure1.names)[1]
    nb_atoms2=size(structure2.names)[1]
    if nb_atoms1 != nb_atoms2
        print("Comparing two structures with different number of atoms. Stopping now!")
        return -1
    end
    rmsd=0
    for atom=1:nb_atoms1
        for i=1:3
            dist=structure1.positions[atom,i]-structure2.positions[atom,i]
            rmsd += dist*dist
        end
    end
    return rmsd
end

end
