include("contactmatrix.jl")

Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]

current_volume=parse(Float64,Volumes[1])
folder=string("/home/moogmt/",current_volume,"/")
file=string(folder,"TRAJEC_wrapped.xyz")

atoms = filexyz.read(file)
