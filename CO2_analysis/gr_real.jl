GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")
push!(LOAD_PATH, GPfolder)

using atom_mod
using cell_mod
using cube_mod
using filexyz
using clustering
using markov
using conversion
using fftw

function computeGr( file_in::T1, V::T2, rmin::T3, rmax::T4, dr::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real }

    # Reading trajectory
    traj, test = filexyz.readFastFile(file_in)
    if ! test
        return zeros(1,1), test
    end

    #Building cell
    cell=cell_mod.Cell_param(V,V,V)

    # Getting nb atoms and steps
    nb_atoms=size(traj[1].names)[1]
    nb_step=size(traj)[1]
    nb_box=Int((rmax-rmin)/dr)
    gr=zeros(nb_box)

    # Computing g_r
    for step=1:nb_step
        print("G(r) Progress: ",step/nb_step*100,"\n")
        for atom1=1:nb_atoms
            for atom2=atom1+1:nb_atoms
                dist=cell_mod.distance(traj[step],cell,atom1,atom2)
                if dist < rmax
                    gr[ Int( trunc( (dist-rmin)/dr ) )+1  ] += 1
                end
            end
        end
    end
    # Normalization
    for i=1:nb_box
        gr[i] = gr[i]*V^3/(nb_step*nb_atoms*(nb_atoms-1)/2)/(4*pi*dr*(i*dr+rmin)^2)
    end

    return gr,test
end

function computeGr( file_in::T1, file_out::T2, V::T3, rmin::T4, rmax::T5, dr::T6 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }

    gr,test=computeGr(file_in,V,rmin,rmax,dr)

    if ! test
        return zeros(1,1),false
    end

    file_o=open(file_out,"w")
    for i=1:size(gr)[1]
        Base.write(file_o,string(i*dr+rmin," ",gr[i],"\n"))
    end
    close(file_o)

    return gr, test
end

function computeFQ( gr::Vector{T1}, rmin::T2, rmax::T3, dr::T4, rho::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real }
    nb_box=1000
    nb_box_int=Int(trunc((rmax-rmin)/dr))
    fq=zeros(nb_box)
    q=zeros(nb_box)
    dq=0.01
    for i=1:nb_box
        loc=0
        q[i] = i*dq
        for j=1:nb_box_int-1
            r=j*dr+rmin
            r2=(j+1)*dr+rmin
            loc += (sin(q[i]*r)*r*(gr[j]-1)+sin(q[i]*r2)*r2*(gr[j+1]-1))*0.5*dr
        end
        loc *= 4*pi/q[i]*rho
        fq[i] = 1 + loc
    end
    return q, fq
end

function computeFQ( file_out::T1, gr::Vector{T2}, rmin::T3, rmax::T4, dr::T5, rho::T6 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    q,fq=computeFQ(gr,rmin,rmax,dr,rho)

    file_o=open(file_out,"w")
    for i=1:size(fq)[1]
        if q[i] > 0
            Base.write(file_o,string(q[i]," ",fq[i],"\n"))
        end
    end
    close(file_o)

    return q,fq
end

# Folder for data
#folder_base="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"
folder_base="/home/moogmt/Data/CO2/CO2_AIMD/"


Temperatures=[2000,3000]
Volumes=[9.8,8.82]

for V in Volumes
    for T in Temperatures

        folder_in=string(folder_base,V,"/",T,"K/")
        file_in=string(folder_in,"TRAJEC_wrapped.xyz")
        folder_out=string(folder_in,"Data/")

        file_out_gr=string(folder_out,"gr.dat")

        rmin=0
        rmax=V/2
        dr=0.001
        gr,test=computeGr( file_in, file_out_gr, V, rmin, rmax, dr )

        file_out_fq=string(folder_out,"fq.dat")

        rho=96/V^3

        q,fq=computeFQ(file_out_fq,gr,rmin,rmax,dr,rho)

    end
end
#
#
# rho=2
#
# n_sk=1000
# sk=ones(n_sk)
# dk=0.05
#
# for i=1:n_sk
#     int_=0
#     for j=1:nb_box
#         int_ += (j*dr+rmin)^2*gr[j]*sin((i*dk*(j*dr+rmin)))/((i*dk*(j*dr+rmin)))*dr
#     end
#     sk[i] += 4*pi*rho*int_
# end
#
# file_out=open(string(folder_out,"sk_real.dat"),"w")
# for i=1:nb_box
#     write(file_out,string(i*dk," ",sk[i],"\n"))
# end
# close(file_out)
#
# y=zeros(20000)
# x=zeros(20000)
# dx=0.05
#
# for i=1:1000
#     x[i] = i*dx
#     y[i] =sin(i*dx)
# end
# ddx=dct(x)
# ddy=dct(y)
#
# file_out=open(string(folder_out,"test_truc.dat"),"w")
# for i=1:nb_box
#     write(file_out,string(x[i]," ",y[i]," ",ddx[i]," ",ddy[i]," \n"))
# end
# close(file_out)
