include("contactmatrix.jl")
include("xyz.jl")

volumes=[8.82,9.0,9.05,9.1,9.2,9.3,9.35,9.4,9.5,9.8]
temperatures=[2000,2500,3000]

V=8.82
T=3000

ps2fs=0.001
timestep=0.5
stride_sim=5
unit=ps2fs*timestep*stride_sim
start_time=5
start_step=Int(start_time/unit)

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder,"TRAJEC.xyz")

print("Reading XYZ file\n")
atoms = filexyz.read(file,stride,start_step)
