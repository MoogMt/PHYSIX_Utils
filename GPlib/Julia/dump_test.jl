include("contactmatrix.jl")

func="PBE-MT"
temperature=3000
volume=[8.82,9.0,9.05,9.1,9.2,9.3,9.35,9.375,9.4,9.5,9.8]
damp_param_list=[20]

T=temperature

for V in volume
    #V=volume[2]

    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/",func,"/",V,"/",T,"K/")
    #folder=string("/home/moogmt/CO2/CO2_AIMD/",V,"/",T,"K/")
    file_in=string(folder,"TRAJEC_wrapped.xyz")

    function sortmatrixindex( x )
        sizex=size(x)[1]
        indexes=zeros(sizex)
        for i=1:sizex
            indexes[i]=i
        end
        for i=1:sizex
            for j=i:sizex
                if x[i] > x[j]
                    stock=x[i]
                    stock2=indexes[i]
                    indexes[i]=indexes[j]
                    x[i]=x[j]
                    indexes[j]=stock2
                    x[j]=stock
                end
            end
        end
        return x,indexes
    end

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
    #damp_param=1
    exchange_param=10

    atoms = filexyz.read( file_in, stride_analysis, start_step )
    cell=cell_mod.Cell_param( V, V, V )

    nb_steps=size(atoms)[1]
    nb_atoms=size(atoms[1].names)[1]

    nb_step_mod=5000
    distance3=zeros(nb_step_mod*nbC)
    angles1=zeros(nb_step_mod*nbC)
    angles2=zeros(nb_step_mod*nbC)
    file=open(string("/home/moogmt/Stupid_test-",V,"-",T,".dat"),"w")
    for step=1:nb_step_mod
        print(string("Progress: ",step/nb_step_mod*100,"%\n"))
        for carbon=1:nbC
            distances=zeros(nbO)
            for oxygen=nbC+1:nbC+nbO
                distances[oxygen-nbC] = cell_mod.distance(atoms[step],cell,carbon,oxygen)
            end
            distances,index=sortmatrixindex(distances)
            distance3[(step-1)*nbC+carbon] = distances[3]
            write(file,string(step*unit," "))
            for i=1:4
                write(file,string(distances[i]," "))
            end
            a=cell_mod.distance(atoms[step],cell,carbon,Int(index[1]+nbC))
            b=cell_mod.distance(atoms[step],cell,carbon,Int(index[2]+nbC))
            c=cell_mod.distance(atoms[step],cell,Int(index[1]+nbC),Int(index[2]+nbC))
            write(file,string( acosd((a*a+b*b-c*c)/(2*a*b))," ") )
            angles1[(step-1)*nbC+carbon] = acosd((a*a+b*b-c*c)/(2*a*b))
            a=cell_mod.distance(atoms[step],cell,carbon,Int(index[2]+nbC))
            b=cell_mod.distance(atoms[step],cell,carbon,Int(index[3]+nbC))
            c=cell_mod.distance(atoms[step],cell,Int(index[2]+nbC),Int(index[3]+nbC))
            write(file,string( acosd((a*a+b*b-c*c)/(2*a*b))," ") )
            write(file,"\n")
            angles2[(step-1)*nbC+carbon] = acosd((a*a+b*b-c*c)/(2*a*b))
        end
    end
    close(file)

    hist2D=zeros(91,91)
    for carbon=1:nb_step_mod*nbC
        print(string("Progress: ",carbon/(nb_step_mod*nbC)*100,"%\n"))
        for angle1=90:180
            for angle2=90:180
                if angles1[carbon] > angle1-0.5 && angles1[carbon] < angle1+0.5
                    if angles2[carbon] > angle1-0.5 && angles2[carbon] < angle1+0.5
                        hist2D[angle1-90+1,angle2-90+1]+=1
                    end
                end
            end
        end
    end

    file=open(string("/home/moogmt/AnglesCO3-",V,"-",T,".dat"),"w")
    for i=1:91
        for i=1:91
            write(file,string(i-1," ",j-1," ",hist2D[i,j],"\n"))
        end
        write(file,"\n")
    end
    close(file)

end
