include("contactmatrix.jl")

function sortmatrix( x )
    sizex=size(x)[1]
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                x[i]=x[j]
                x[j]=stock
            end
        end
    end
    return x
end

temperature=[2000,2250,25000,3000]
volume=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.325,9.35,9.375,9.4,9.5,9.8]

# Volume loop
for T in temperatures
for V in volume
    #V=9.1
    # Temp
    #T=temperature[1]

    # Trajectory file
    folder=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/",V,"/",T,"K/")
    file_in=string(folder,"TRAJEC_wrapped.xyz")

    if ! isfile(file_in)
        print("\n")
        break
    end

    # Time units
    stride_sim=5
    fs2ps=0.001
    time_sim=1 # in fs
    unit=time_sim*fs2ps*stride_sim# in ps

    # Stride and step
    stride_analysis=1
    start_time=5
    start_step=Int(start_time/(unit*stride_sim))

    # Reading trajectory and cell
    atoms = filexyz.read( file_in, stride_analysis, start_step )
    cell=cell_mod.Cell_param( V, V, V )

    # Number of aotms and steps
    nb_steps=size(atoms)[1]
    nb_atoms=size(atoms[1].names)[1]

    # Computing NN
    nn=zeros(nb_steps*32,5)
    for step=1:nb_steps
        print("progress: ",step/nb_steps*100,"%\n")
        for carbon=1:32
            distances=zeros(nb_atoms-32)
            for oxygen=33:nb_atoms
                distances[oxygen-32] = cell_mod.distance(atoms[step],cell,carbon,oxygen)
            end
            distances=sortmatrix(distances)
            nn[(step-1)*32+carbon,1]=step
            for i=2:5
                nn[(step-1)*32+carbon,i]=distances[i-1]
            end
        end
    end

    dfrac=0.05
    folder="/home/moogmt/"
    time_window=0.5 # in ps
    for k=3:4

        hist_bond=zeros(Int(trunc(nb_steps/(time_window/unit))))
        hist_unbond=zeros(Int(trunc(nb_steps/(time_window/unit))))
        hist_trans=zeros(Int(trunc(nb_steps/(time_window/unit))))

        nb_steps=Int(trunc(nb_steps/(time_window/unit)))*(time_window/unit)

        count_step=0
        count = 1

        file_out=open(string(folder,"Signal-",k,"_V-",V,"_T-",T,".dat"),"w")
        for step=1:nb_steps
            count_bond=0
            count_unbond=0
            for carbon=1:32
                if nn[Int((step-1)*32+carbon),k+1] < 1.5 && nn[Int((step-1)*32+carbon),k+1]  > 1.3
                    count_bond += 1
                elseif nn[Int((step-1)*32+carbon),k+1] > 1.9 && nn[Int((step-1)*32+carbon),k+1] < 2.3
                    count_unbond += 1
                end
            end
            if count_step == Int(time_window/unit)
                norm=(count_bond+count_unbond)
                string2=string(count*(time_window/unit)*unit," ",count_bond/norm," ",count_unbond/norm,"\n")
                write(file_out,string2)
                hist_bond[count] = count_bond/norm
                hist_unbond[count] = count_unbond/norm
                count += 1
                count_step = 0
            else
                count_step += 1
            end

        end
        close(file_out)

        # Hist2D
        if k == 3

            # Aggregating data
            hist2D=zeros(Int(1/dfrac),Int(1/dfrac))
            for i=1:1:nb_steps/(time_window/unit)
                check=false
                for j=1:1/dfrac
                    for k=1:1/dfrac
                        if hist_bond[Int(i)] < j*dfrac  &&  hist_bond[Int(i)] > (j-1)*dfrac && hist_unbond[Int(i)] < k*dfrac && hist_unbond[Int(i)] > (k-1)*dfrac
                            hist2D[Int(j),Int(k)] += 1
                            check = true
                            break
                        end
                    end
                    if check == true
                        break
                    end
                end
            end

            # Writting data
            file_out2=open(string(folder,"HistSign-V-",V,"_T-",T,".dat"),"w")
            for i=1:1/dfrac
                for j=1:1/dfrac
                    string2=string(i*dfrac," ",j*dfrac," ",hist2D[Int(i),Int(j)],"\n")
                    write(file_out2,string2)
                end
                write(file_out2,"\n")
            end
            close(file_out2)

        end
        # End of Hist2D

    end
    # End of k
end
# End of V
end
# End of T
