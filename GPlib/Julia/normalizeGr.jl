# Loading file
include("contactmatrix.jl")


Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[1750,2000,2250,2500,2750,3000]

folder_base=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

for V in Volumes
    for T in Temperatures
        folder_local = string( folder_base , V , "/" , T , "K/Data/" )
        if isfile( string(folder_local,"gCO.dat") )
            file_out = open( string( folder_local, "gCO_norm.dat" ), "w" )
            # Read
            #--------------------------------------------
            file  = open( string( folder_local, "gCO.dat" ) )
            lines = readlines( file )
            close( file )
            #--------------------------------------------
            r  = zeros( size( lines )[1] )
            gr = zeros( size( lines )[1] )
            max=0
            for i=1:size( lines )[1]
                r[i]  = parse( Float64, split( lines[i] )[1] )
                gr[i] = parse( Float64, split( lines[i] )[2] )
                if gr[i] > max
                    max=gr[i]
                end
            end
            gr /= max
            for i=1:size(lines)[1]
                write( file_out , string( r[i] , " " , gr[i] , "\n" ) )
            end
            close( file_out )
        end
    end
end
