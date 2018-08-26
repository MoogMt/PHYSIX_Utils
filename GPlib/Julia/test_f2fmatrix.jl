points=zeros(4,2)

points[2,1]=1
points[3,2]=1
points[4,:]=[1,1]

folder="/home/moogmt/CO2/CO2_AIMD/clusters/"

dmax=0
dist=zeros(4,4)
for i=1:3
    for j=i+1:4
        dist[i,j]=sqrt( sum( (points[i,:]-points[j,:]).*(points[i,:]-points[j,:]) ) )
        dist[j,i]=dist[i,j]
        if dmax < dist[i,j]
            dmax = dist[i,j]
        end
    end
end

dist=dist/dmax

data=open(string(folder,"FRAME_TO_FRAME.MATRIX"),"w")
write(data,string(4," ",dmax,"\n"))
for i=1:4
    for j=1:4
        write(data,string(dist[i,j]," "))
    end
    write(data,string("\n"))
end
close(data)
