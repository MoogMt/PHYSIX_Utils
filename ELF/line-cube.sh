awk 'BEGIN{ ix=1; iy=1; iz=1; ddcut=2.*2. }
{
 if (FILENAME=="list") {
   nc=$NF
   for (ic=1;ic<=nc;ic++) { iat[ic]=$ic }
 }else{
   nn++
   #if (nn<=2||nn<=6+nat) print
   if (nn==3) {
     nat=$1
     for (k=1;k<=3;k++) { orig[k]=$(k+1); }
   }
   if (nn>=4&&nn<=6) {
     i++
     nvox[i]=$1
     for (k=1;k<=3;k++) { delta[i,k]=$(k+1) }
     for (ix=0;ix<nvox[1];ix++) {
       for (iy=0;iy<nvox[2];iy++) {
         for (iz=0;iz<nvox[3];iz++) {
           ivox++
           # it works only for orthorombic cube vectors
           posvox[ivox,1]=orig[1]+ix*delta[1,1]
           posvox[ivox,2]=orig[2]+iy*delta[2,2]
           posvox[ivox,3]=orig[3]+iz*delta[3,3]
         }
       }
     }
     ivox=0
   }
   if (nn>=7&&nn<=6+nat) {
     j++
     for (k=1;k<=3;k++) { pos[j,k]=$(k+2) }
   }

   if (nn==6+nat) {
#     print "**********************"
     nv=20
     for(i=1;i<=nv;i++){
       for (h=1;h<=3;h++) { v[i,h]=pos[iat[1],h]+(pos[iat[2],h]-pos[iat[1],h])*i/nv
       }
#         print v[i,1],v[i,2],v[i,3]
     }
#     print pos[iat[1],h],pos[iat[2],h]
     ddmax=(delta[1,1]+delta[2,2]+delta[3,3])/3.
     ddmax*=ddmax
     print ddmax
#     print "**********************"
   }

   if (nn>6+nat) {
     for (k=1;k<=NF;k++) {
       ivox++
       rho=$k
       ok=0
#       for (ic=1;ic<=nc;ic++) {
#         dd=0.
#         for (h=1;h<=3;h++) {
#           # warning: no pbc are used here
#           dx=pos[iat[ic],h]-posvox[ivox,h]
#           dd+=dx*dx
#         }
#         if (dd<=ddcut) ok=1
#       }
#       if (ok==0) rho=0.
#       printf "%13.5e",rho

       for(i=1;i<=nv;i++){
         dd=0.
         for (h=1;h<=3;h++) {
           dd+=(v[i,h]-posvox[ivox,h])*(v[i,h]-posvox[ivox,h])
         }
         if (dd<ddmax) print i,rho > "ELF_on_line"
         #print i,dd,rho > "ELF_on_line"
       }

     }
     #printf "\n"
   }
 }
}' list $1
