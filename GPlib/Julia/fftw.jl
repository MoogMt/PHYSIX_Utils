module fftw

using FFTW

function doFourierTransform( signal::Vector{T1}, dt ) where { T1 <: Real }
    return fftfreq(size(signal)[1],1/dt), abs.(fft(signal) )
end

function doFourierTransformShift( signal::Vector{T1}, dt ) where { T1 <: Real }
    return fftfreq(size(signal)[1],1/dt) |> fftshift, abs.(fft(signal) |> fftshift )
end

end

#  Examples
# dt=0.001
# x=range(0,stop=5*pi,step=dt)
# y=cos.(2*pi*x*50)
#
# j,k=doFourierTransform(y,dt)
#
# file_out=open(string("/home/moogmt/test2.dat"),"w")
# for i=1:size(z)[1]
#     Base.write(file_out,string(x[i]," ",y[i],"\n"))
# end
# close(file_out)
#
# file_out=open(string("/home/moogmt/test.dat"),"w")
# for i=1:size(k)[1]
#     Base.write(file_out,string(j[i]," ",k[i],"\n"))
# end
# close(file_out)
