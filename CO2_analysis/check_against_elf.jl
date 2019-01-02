GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))
include(string(GPfolder,"cubefile.jl"))
include(string(CO2folder,"markovCO2.jl"))

folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/ELF/8.82/Trajectory_2/"

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75
cut_off_states = 0.1 # 0.1% cut-off for a state to be considered statististically viable

min_lag=1    # min tau
max_lag=5001 # max tau
d_lag=5      # delta tau
unit=0.005   # units of the simulation

V=8.82
nb_steps=100

atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,1,"_elf.cube") )
cell=cell_mod.Cell_param(cell_mod.cellMatrix2Params(cell_matrix))
atoms=atom_mod.move(atoms,(-1)*elf.origin )
atoms=cell_mod.wrap(atoms,cell)
elf_value=cube_mod.dataInTheMiddleWME( atoms, cell , 25, 35, elf )

print("Computing Data\n")
for i=1:nb_steps
    atoms, cell_matrix, elf = cube_mod.readCube( string(folder_base,i,"_elf.cube") )
    density = cube_mod.readCube( string(folder_base,i,"_elf.cube") )[3]

end

data=buildCoordinationMatrix( traj , cell , cut_off_bond )


state_matrix, percent, unused_percent = assignDataToStates( data , states )

nb_states=size(states)[1]

coordinances=zeros(nb_states)
for i=1:nb_states
    for j=1:nb_dim
        if states[i,j] > 0
            coordinances[i] += 1
        end
    end
end
