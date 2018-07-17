include("contactmatrix.jl")
include("pdb.jl")

func="PBE-MT"
temperature=3000
volume=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8]


V=volume[1]
T=temperature

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
file_in=string(folder,"TRAJEC_wrapped.xyz")

# Sim parameters
stride_sim=5
fs2ps=0.001
time_sim=0.5 # in fs
unit=time_sim*fs2ps*stride_sim# in ps
stride_analysis=1
start_time=5
start_step=Int(start_time/(unit*stride_sim))
nbC=32
nbO=64
cut_off=1.6

# Reading trajectory
atoms = filexyz.read( file_in, stride_analysis, start_step )
cell=cell_mod.Cell_param( V, V, V )

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

atoms=atoms[1:10:nb_steps]

nb_steps=size(atoms)[1]

file=open(string("/home/moogmt/trajec-",V,"-",T,"-"func,".pdb"),"w")
for i=1:nb_steps
    pdb.writeStep(atoms[i],cell,file)
end
close(file)
