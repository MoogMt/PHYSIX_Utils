module cube_mod

include("atom.jl");
include("cell.jl")

# Reads a cube file and returns all or parts of its informations
function readCube{}()

    
    return nb_atoms, atom_list, density
end

end
