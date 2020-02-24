awk 'BEGIN{ niter=1000; bestmaxd=0.; dmax=0 }
{
 if (NR==1) print "reading map ..."
 if(NF>0 && FILENAME=="map"){ d[$1,$2]=$3; if($3>dmax)dmax=$3; n=$1 }
# if(FILENAME=="restart"){ 
#  m++
#  x[m]=$1; y[m]=$2
# }
}END{
 print "n =",n
 #n=13
 k=1.

 srand()

 for(i=1;i<=n;i++){
  x[i]=rand()*dmax
  y[i]=rand()*dmax
 }
 
 for(it=1;it<=niter;it++){

   for(i=1;i<=n;i++){
     sx[i]=0
     sy[i]=0
   }

   e=0.
   maxd=0.
   for(i=1;i<n;i++){
    for(j=i+1;j<=n;j++){
     dxij=x[i]-x[j]
     dyij=y[i]-y[j]
     dij=dxij*dxij+dyij*dyij
     if (dij==0) dij=0.000001
     dij=sqrt(dij)
     if (it==niter) printf "%8.5f %8.5f\n",d[i,j],dij > "errors"
     delta=dij-d[i,j]
     delta2=delta*delta
     if(delta2>maxd){ maxd=delta2 }
     e+=0.5*k*delta2
     sx[i]+=-k*delta*dxij/dij
     sx[j]-=-k*delta*dxij/dij
     sy[i]+=-k*delta*dyij/dij
     sy[j]-=-k*delta*dyij/dij
    } 
   }

   if(it%1==0){
     printf "progress = %6.4f   energy = %10.7f   maxdev = %8.5f\n",it/niter,e,sqrt(maxd)
   }

   smax=0
   for(i=1;i<=n;i++){
     if (sx[i]>smax) smax=sx[i]
     if (sy[i]>smax) smax=sy[i]
   }
   scale=1.0
   if (smax>dmax*0.05) scale=dmax*0.05/smax
   for(i=1;i<=n;i++){
     x[i]+=sx[i]*scale
     y[i]+=sy[i]*scale
   }

 }

 printf "\n"
 print "final energy =",e
 print "average and max deviation from real distances =",sqrt(2.*e/(k*n*(n-1)/2.)),sqrt(maxd)

 print "writing landscape points in file points"
 for(i=1;i<=n;i++){ printf "%8.3f %8.3f %4d\n",x[i],y[i],i > "points" }

}' map 
