include("contactmatrix.jl")

Volumes=["8.82","9.0","9.05","9.1","9.15","9.2","9.3","9.35","9.4","9.8"]
Temperature=[2000,2250,2500,3000,3500]

current_volume=parse(Float64,Volumes[4])
current_temperature=Temperature[4]
folder=string("/home/moogmt/",current_volume,"/",current_temperature,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

atoms = filexyz.read(file,20)
