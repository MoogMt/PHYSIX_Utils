GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
CO2folder=string("/home/moogmt/PHYSIX_Utils/CO2_analysis/")

include(string(CO2folder,"markovCO2.jl"))
include(string(GPfolder,"cell.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Number of atoms
nbC=32
nbO=nbC*2
nb_atoms=96
V=8.82
T=3000
cell=cell_mod.Cell_param(V,V,V)

folder_out=string(folder_base,V,"/",T,"K/")

file_in=open(string(folder_base,V,"/",T,"K/1-run/FTRAJECTORY2"))
lines=readlines(file_in)
close(file_in)
nb_steps=Int(trunc(size(lines)[1]/96))
positions=zeros(Real,nb_steps,nb_atoms,3)
forces=zeros(Real,nb_steps,nb_atoms,3)
for i=1:nb_steps
    for j=1:nb_atoms
        line=split(lines[(i-1)*nb_atoms+j])
        for k=1:3
            positions[i,j,k] = cell_mod.wrap(parse( Float64, line[k+1] ),cell.length[k])
        end
    end
end

file_energy=open(string(folder_base,V,"/",T,"K/1-run/ENERGIES"))
lines=readlines(file_energy)
close(file_energy)

nb_energy=size(lines)[1]
energy=zeros(nb_energy)
for i=1:nb_energy
    energy[i]=parse(Float64,split(lines[i])[3])
end
energy=energy[1:5:nb_energy]


cell_matrix=zeros(3,3)
for i=1:3
    cell_matrix[i,i]=V
end

atom_names=Vector{AbstractString}(undef,nb_atoms)
for i=1:nbC
    atom_names[i]="C"
end
for i=33:nbC+nbO
    atom_names[i]="O"
end

charge=zeros(Real,nb_atoms)
for i=1:nbC
    charge[i]=6
end
for i=33:nbC+nbO
    charge[i]=8
end

charge_system=zeros(nb_steps)

function writeN2P2conf( positions::Array{T1,3}, forces::Array{T2,3}, cell_matrix::Array{T3,2}, atom_types::Vector{T4}, charge::Vector{T5}, charge_system::Vector{T6}, energy::Vector{T7}, file_out::T8 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: AbstractString, T5 <: Real, T6 <: Real, T7<:Real, T8 <: AbstractString }
    nb_steps,nb_atoms,dim=size(positions)
    file=open(file_out,"w")
    for step=1:nb_steps
        write(file,string("begin\n"))
        write(file,string("comment STEP ",step,"\n"))
        for i=1:3
            write(file,string("lattice "))
            for j=1:3
                write(file,string(cell_matrix[i,j]," "))
            end
            write(file,string("\n"))
        end
        for atom=1:nb_atoms
            write(file,string("atom "))
            for i=1:3
                write(file,string(positions[step,atom,i]," "))
            end
            write(file,string(atom_types[atom]," "))
            write(file,string(charge[atom]," "))
            write(file,string("X "))
            for i=1:3
                write(file,string(forces[step,atom,i]," "))
            end
            write(file,string("\n"))
        end
        write(file,string("energy ",energy[step],"\n"))
        write(file,string("charge ",charge_system[step],"\n"))
        write(file,string("end\n"))
    end
    close(file)
end

writeN2P2conf(positions,forces,cell_matrix,atom_names,charge,charge_system,energy,"/home/moogmt/input.data")
