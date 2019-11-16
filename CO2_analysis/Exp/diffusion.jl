GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using clustering
using filexyz
using pdb
using markov
using fftw
using correlation
using conversion
using exp_data
volumes=[9.1,9.2,9.3,9.4]
#temperatures=[2000,2500,3000]
temperatures=[3000]

# for V in volumes
# for V in volumes
# for T in temperatures
V=8.82
T=3000

unit=0.0005*5
stride=1
start_time=5
start_step=Int(start_time/unit)

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder,"TRAJEC.xyz")

print("Reading XYZ file\n")
atoms = filexyz.read(file,stride,start_step)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

# File
output=string("/home/moogmt/diff_",V,"_",T,"K_.dat")
file=open(output,"w")

C_MSD=zeros(nb_steps)
O_MSD=zeros(nb_steps)
All_MSD=zeros(nb_steps)

for step=1:nb_steps
    print("V = ",V," , T = ",T," , progress: ",step/nb_steps*100,"%\n")

    #---------------------------------------
    # Determination of gravicenter
    #---------------------------------------
    local_gc=zeros(3)
    for index=1:nb_atoms
        local_gc += atoms[step].positions[index,:]
    end
    local_gc /= nb_atoms
    #---------------------------------------

    #---------------------------------------
    # Bookeep initial GC
    #---------------------------------------
    for index=1:nb_atoms
        atoms[step].positions[index,:] -= local_gc
    end
    #---------------------------------------
    # Compute MSD
    #---------------------------------------
    for index=1:32
        tot=0
        for i=1:3
            tot += (atoms[step].positions[index,i]-atoms[1].positions[index,i])^2
        end
        C_MSD[step] += tot
    end
    for index=33:nb_atoms
        tot=0
        for i=1:3
            tot += (atoms[step].positions[index,i]-atoms[1].positions[index,i])^2
        end
        O_MSD[step] += tot
    end
    #---------------------------------------

    #---------------------------------------
    # Normalizing
    #---------------------------------------
    All_MSD[step] = (C_MSD[step]+O_MSD[step])/nb_atoms
    C_MSD[step] /=32
    O_MSD[step] /= 64
    #---------------------------------------

    #---------------------------------------
    # Writting results
    #---------------------------------------
    out_string=string(step*unit*stride," ",C_MSD[step]," ",O_MSD[step]," ",All_MSD[step],"\n")
    write(file,out_string)
    #---------------------------------------
end
close(file)

file_check=open(string("/home/moogmt/varMSD_",V,"_",T,".dat"),"w")
for i=2:nb_steps
    dC=abs((C_MSD[i]-C_MSD[i-1])/(unit*stride))
    dO=abs((O_MSD[i]-O_MSD[i-1])/(unit*stride))
    write(file_check,string(i*unit*stride," ",dO-dC,"\n"))
end
close(file_check)


# end
# end
# end
