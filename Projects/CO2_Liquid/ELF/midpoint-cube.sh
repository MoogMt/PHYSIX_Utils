awk 'BEGIN{ ix=1; iy=1; iz=1; ddcut=2.*2. } 
{
 if (FILENAME=="list") {
   nc=$NF
   for (ic=1;ic<=nc;ic++) { iat[ic]=$ic }
 }else{
   nn++
   if (nn<=2||nn<=6+nat) 
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
     for (h=1;h<=3;h++) 
     { 
        v[h]=(pos[iat[2],h]+pos[iat[1],h])/2 
     }
     ddmax=(delta[1,1]+delta[2,2]+delta[3,3])/3.
     ddmax*=ddmax
    }

   if (nn>6+nat) {
     for (k=1;k<=NF;k++) { 
       ivox++
       rho=$k
       dd=0.
       for (h=1;h<=3;h++) {
          dd+=(v[h]-posvox[ivox,h])*(v[h]-posvox[ivox,h])
       }
       if (dd<ddmax) print rho > "elf_midpoint"
     }
   }
 }
}' list $1
