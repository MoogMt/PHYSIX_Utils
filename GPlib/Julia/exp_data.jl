module exp_data

using atom_mod
using cell_mod
using filexyz
using pdb
using fftw
using correlation
using conversion

# Vibration Density of States
#==============================================================================#
function vdosFromPosition( file_traj::T1 , max_lag_frac::T2 , to_nm::T3, dt::T4 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real }

        # Reading Trajectory
        traj,test=readFastFile(file_traj)

        if ! test
            return zeros(1,1), zeros(1,1), test
        end

        # Computing velocities
        velocity=cell_mod.velocityFromPosition(traj,dt,dx)

        nb_atoms=size(velocity)[2]
        nb_step=size(velocity)[1]

        # Compute scalar product
        velo_scal=zeros(nb_step,nb_atoms)
        for atom=1:nb_atoms
            for step=1:nb_step
                for i=1:3
                    velo_scal[step,atom] += velocity[step,atom,i]*velocity[1,atom,i]
                end
            end
        end

        max_lag=Int(trunc(nb_step*max_lag_frac))

        # Average correlation
        freq=zeros(max_lag)
        vdos=zeros(max_lag)
        for atom=1:nb_atoms
            freq,vdos_loc=fftw.doFourierTransformShift( correlation.autocorrNorm( velo_scal[:,atom] , max_lag ), dt )
            vdos += vdos_loc
        end
        vdos /= nb_atoms

        # Conversion to cm-1
        freq=freq.*conversion.Hz2cm

    return freq, vdos, test
end

function vdosFromPosition( file_traj::T1 , file_out::T2 , max_lag_frac::T3 , to_nm::T4, dt::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real, T5 <: Real }

    freq,vdos,test=vdosFromPosition( file_traj , max_lag_frac , to_nm, dt )
    if ! test
        return zeros(1,1),zeros(1,1),false
    end

    # Writting data to file
    file_o=open(string(file_out),"w")
    for i=1:size(vdos)[1]
        if freq[i] > 0 # Remove frequency 0 and symmetric
            Base.write(file_o,string(freq[i]," ",vdos[i],"\n"))
        end
    end
    close(file_o)

    return freq, vdos, test
end
#==============================================================================#

# Radial Distribution Function
#==============================================================================#
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
#==============================================================================#

# Structure Factor
#==============================================================================#
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
#==============================================================================#

end
