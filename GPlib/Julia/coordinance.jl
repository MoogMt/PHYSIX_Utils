include("contactmatrix.jl")
include("xyz.jl")

volumes=[9.1,9.2,9.3,9.4]
volumes=[8.82,9.0,9.8]
temperatures=[3000]

unit=0.0005
unit=0.0005*5
stride=5
stride=1
# Cutting the first 5 ps
start_time=5
start_step=Int(start_time/unit)
cut_off=1.6

#V=volumes[1]
for V in volumes
T=temperatures[1]

folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
file=string(folder,"TRAJEC_wrapped.xyz")

cell=cell_mod.Cell_param(V,V,V)

atoms = filexyz.read(file,stride,start_step)

nb_steps=size(atoms)[1]
nb_atoms=size(atoms[1].names)[1]

coord2=zeros(nb_steps)
coord3=zeros(nb_steps)
coord4=zeros(nb_steps)


file_coord=open(string("/home/moogmt/coord_frac_",T,"_",V,"_",cut_off,".dat"),"w")

for step=1:nb_steps
    print("progress: ",step/nb_steps*100,"% \n")
    coord_avg=0
    coord2=0
    coord3=0
    coord4=0
    for atom=1:32
        coord=0
        for atom2=33:nb_atoms
            if atom != atom2
                if cell_mod.distance(atoms[step],cell,atom,atom2) < cut_off
                    coord += 1
                end
            end
        end
        if coord == 4
            coord4 += 1
        elseif coord == 3
            coord3 += 1
        elseif coord == 2
            coord2 += 1
        end
        coord_avg += coord
    end
    write(file_coord,string(step*unit*stride," ",coord2/32," ",coord3/32," ",coord4/32," ",coord_avg/32,"\n"))
end
close(file_coord)
end


#==============================================================================#

volumes   = [ 8.82, 9.0, 9.1, 9.2, 9.3, 9.4, 9.8 ]
pressures = [ 65,   60,   55,  53,  50,  45,  30 ]
T=3000

cut_offs=[1.6,1.7,1.85]

for cutoff in cut_offs

file_coord3=open(string("/home/moogmt/coord_avg",cutoff,".dat"),"w")
count=1

for V in volumes

file_coord2=open(string("/home/moogmt/coord_frac_",T,"_",V,"_",cutoff,".dat"))
lines=readlines(file_coord2)
close(file_coord2)

nb_points=size(lines)[1]

avg2=0
avg3=0
avg4=0
for i=1:nb_points
    line=split(lines[i])
    avg2 += parse( Float64, line[2] )
    avg3 += parse( Float64, line[3] )
    avg4 += parse( Float64, line[4] )
end
avg2 /= nb_points
avg3 /= nb_points
avg4 /= nb_points
write(file_coord3,string(pressures[count]," ",avg2," ",avg3," ",avg4,"\n"))
count+=1
end

close(file_coord3)

end
