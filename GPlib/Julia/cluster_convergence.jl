# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=10
co=cut_off[1]

nbC=32
nbO=2*nbC
