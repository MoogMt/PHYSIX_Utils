include("contactmatrix.jl")
include("xyz.jl")

volumes=[8.82,9.0,9.05,9.1,9.2,9.3,9.4,9.8]
temperatures=[2000,2500,3000]

function sortmatrix( x )
    sizex=size(x)[1]
    index_x=zeros(sizex)
    for i=1:sizex
        index_x[i]=i
    end
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stocki=index_x[i]
                x[i]=x[j]
                index_x[i]=index_x[j]
                x[j]=stock
                index_x[j]=stocki
            end
        end
    end
    return index_x
end

T=300#temperatures[1]

#for V in volumes
    V=volumes[1]
    #     for T in temperatures

    ps2fs=0.001
    timestep=0.5
    stride = 1
    unit=ps2fs*timestep*stride
    start_time=5
    start_step=Int(start_time/unit)
    nbC=32
    nbO=2*nbC

    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
    file=string(folder,"TRAJEC_wrapped.xyz")

    if isfile( string(file) )

        print("Reading XYZ file\n")
        traj = filexyz.read(file,stride,start_step)
        cell=cell_mod.Cell_param(V,V,V)

        nb_steps=size(traj)[1]
        nb_atoms=size(traj[1].names)[1]

        fileC=open(string(folder,"/allInfoCvO-",V,"-",T,"K.dat"),"w")
        for step=1:5:nb_steps
            print("Progres: ",step/nb_steps*100,"%\n")
            for carbon=1:nbC
                distances = zeros(nbO)
                for oxygen=1:nbO
                    distances[oxygen] = cell_mod.distance(traj[step],cell,carbon,oxygen+nbC)
                end
                index=sortmatrix( distances )
                write(fileC, string(step*unit," ") )
                for i=1:4
                    write(fileC, string(distances[i]," ") )
                end
                maxN=4
                for i=1:maxN-1
                    for j=i+1:maxN
                        a=cell_mod.distance(traj[step],cell,carbon,Int(index[i]+nbC))
                        b=cell_mod.distance(traj[step],cell,carbon,Int(index[j]+nbC))
                        c=cell_mod.distance(traj[step],cell,Int(index[i]+nbC),Int(index[j]+nbC))
                        write(fileC,string(acosd((a*a+b*b-c*c)/(2*a*b))," "))
                    end
                end
                write(fileC,string("\n"))
            end
        end
        close(fileC)
    else
        print("File ", string(folder,file), " does exists...\n")
    end

    #     end
#end
