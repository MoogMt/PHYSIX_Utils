cat > topol.top << ***
; Include forcefield parameters
;#include "oplsaa.ff/forcefield.itp"

[ defaults ]
; nbfunc (1=LJ)   comb-rule   gen-pairs   fudgeLJ fudgeQQ
1                 1           no          1        1 

[ atomtypes ]
; note: C6  = 4*eps*sigma^6  [kJ/mol*nm^6]
;       C12 = 4*eps*sigma^12 [kJ/mol*nm^12]
; name  at.num  mass     charge   ptype(A=atom)   V(C6)         W(C12)
  Nh    7       14.0067  -1.02        A          4.552576e-3   7.284736e-6 
  Hn    1        1.0080   0.34        A          0.            0.
  Oh    8       15.9994  -1.04        A          2.552e-3      2.5104e-6
  Ho    1        1.0080   0.52        A          0.            0.

[ moleculetype ]
; molname   nrexcl
NH3          3

[ atoms ]
; note: type = atomtypes, resname = molname, atname = (name in geometry file)
;   nr       type   resnr resname  atname   cgnr     charge    mass  
     1       Nh      1    NH3      NH1       1       -1.02    14.00670
     2       Hn      1    NH3      HN1       1        0.34     1.008
     3       Hn      1    NH3      HN2       1        0.34     1.008
     4       Hn      1    NH3      HN3       1        0.34     1.008

[ angles ]
;  i    j    k  func       th0       cth
   2    1    3    1      109.5    292.880   ; taken from OPLS-AA LYS
   2    1    4    1      109.5    292.880   ; taken from OPLS-AA LYS
   3    1    4    1      109.5    292.880   ; taken from OPLS-AA LYS

[ constraints ]
;  at1  at2   funct   dist(nm)
   1    2     2       0.106
   1    3     2       0.106
   1    4     2       0.106

 
[ moleculetype ]
; molname   nrexcl
H2O           3

[ atoms ]
; note: type = atomtypes, resname = molname, atname = (name in geometry file)
;   nr       type   resnr resname  atname   cgnr     charge    mass  
     1       Oh     1     H2O       OH1       1       -1.04    15.99940    
     2       Ho     1     H2O       HO1       1        0.52     1.008
     3       Ho     1     H2O       HO2       1        0.52     1.008

[ angles ]
;  i    j    k  func       th0    cth 
   2    1    3    1      104.52    836.800    ; taken from OPLS-AA LYS

[ constraints ]
;  at1  at2   funct   dist(nm)
   1    2     2       0.096
   1    3     2       0.096

[ system ]
mysystem

[ molecules ]
***

awk '{
 if (NR==1){ a=$1/10; b=$2/10; c=$3/10 }
 else{
  n++
  s[n]=$1; x[n]=$2; y[n]=$3; z[n]=$4
 }
 mx=2
 my=2
 mz=4
}END{
 print "NH3H2O"
 printf "%5d\n",n*mx*my*mz
 iii=0
 res=0
 for(i=0;i<mx;i++){
 for(j=0;j<my;j++){
 for(k=0;k<mz;k++){
  for(h=1;h<=n;h++){
   iii++
   if (s[h]=="NH1") { print "NH3  1" >> "topol.top"; rrr++; resname="NH3" }
   if (s[h]=="OH1") { print "H2O  1" >> "topol.top"; rrr++; resname="H2O"  }
   
   printf "%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n",rrr,resname,s[h],iii,(x[h]+i)*a,(y[h]+j)*b,(z[h]+k)*c
  }
 }}}
 printf "%10.5f%10.5f%10.5f\n",a*mx,b*my,c*mz
}' cryst-coord > start.gro

