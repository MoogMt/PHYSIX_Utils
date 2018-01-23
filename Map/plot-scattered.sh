echo "check"

awk '{
 if(NR==1) { n=$1; dmax=$2 }else {
   for (i=1;i<=n;i++) { print NR-1,i,$i*dmax }
 }
}' FRAME_TO_FRAME.MATRIX > map

echo "check"

awk 'BEGIN{ niter=2000000; bestmaxd=0. }
{
 if(NF>0 && FILENAME=="map"){ d[$1,$2]=$3 }
# if(FILENAME=="restart"){ 
#  m++
#  x[m]=$1; y[m]=$2
# }
}END{

 n=30
 k=0.02
 s=0.005
 kT=0.0000002
 
 srand()

 for(i=1;i<=n;i++){
  x[i]=rand()
  y[i]=rand()
 }
 e=0.
 for(i=1;i<n;i++){
  for(j=i+1;j<=n;j++){
   dij=(x[i]-x[j])**2+(y[i]-y[j])**2
   if(dij>0.){ dij=sqrt(dij) }
   e+=0.5*k*(dij-d[i,j])**2
   printf "error %2d %2d = %8.5f relative = %8.5f orig_d mapped_d = %8.5f %8.5f\n",i,j,dij-d[i,j],(dij-d[i,j])/d[i,j],d[i,j],dij
  } 
 }
 print "initial energy =",e

 for(i=1;i<=n;i++){ xx[i]=x[i]; yy[i]=y[i] }

 for(it=1;it<=niter;it++){

  ii=int(rand()*n)+1
#  for(i=1;i<=n;i++){
   xx[ii]=x[ii]+(rand()-0.5)*s
   yy[ii]=y[ii]+(rand()-0.5)*s
#   print x[ii],xx[ii],y[ii],yy[ii]
#  }
#   i=int(rand()*n)+1
# #  for(i=1;i<=n;i++){
#    xx[i]=x[i]+(rand()-0.5)*s
#    yy[i]=y[i]+(rand()-0.5)*s
# #  }

  ee=0.
  maxd=0.
  for(i=1;i<n;i++){
   for(j=i+1;j<=n;j++){
    dij=(xx[i]-xx[j])**2+(yy[i]-yy[j])**2
    if(dij>0.){ dij=sqrt(dij) }
    dd=(dij-d[i,j])**2
    if(dd>maxd){ maxd=dd }
    ee+=0.5*k*dd
    if(it==niter) { 
      printf "error %2d %2d = %8.5f relative = %8.5f orig_d mapped_d = %8.5f %8.5f\n",i,j,dij-d[i,j],(dij-d[i,j])/d[i,j],d[i,j],dij
      printf "%8.5f %8.5f\n",d[i,j],dij > "errors"
    }
   } 
  }
  
  r=rand()
  de=(ee-e)/kT
#  print de
  if(de>0.){ if(de<50.){ m=exp(-de) }else{ m=0. } } else { m=2.0 }
  nch++
  if (r<m) {
   for(i=1;i<=n;i++){
    x[i]=xx[i]
    y[i]=yy[i]
    e=ee
#    printf "a"
   }
   nach++
   bestmaxd=maxd
  }else{
#   printf "r" 
  }
  if(nch==1000){
   printf "progress = %6.4f   accepted = %6.4f   energy = %10.7f   maxdev = %8.5f\n",it/niter,nach/nch,e,sqrt(bestmaxd)
   nch=0
   nach=0
  }
 }
 printf "\n"
 print "final energy =",e
 print "average and max deviation from real distances =",sqrt(2.*e/(k*n*(n-1)/2.)),sqrt(bestmaxd)

 for(i=1;i<=n;i++){ printf "set label \"S%d\"   at %8.3f , %8.3f font \"arial,16\"\n",i,x[i],y[i] > "gp-points1" }
#print  "sx=0.005; sy=0.005; unset key; p \"<head -n9 gp-points1\" u ($5-s):($7-s) pt 7 ps 2 lt 3, \"<tail -n17 gp-points1\" u ($5-s):($7-s) pt 7 ps 1 lt 1"    > "gp-points2"
}' map 

