awk '{
 
  if (NR==1){ a=$1/10; b=$2/10; c=$3/10 }
  else{
    n++
    s[n]=$1; x[n]=$2; y[n]=$3; z[n]=$4
  }
  }END{
    print "2NH3H2O"
    print "470"
    iii=0
    res=0
    for(h=1;h<=n;h++){
      iii++
      if (s[h]=="NH1") { rrr++; resname="NH3" }
      if (s[h]=="OH1") { rrr++; resname="H2O"  }
        printf "%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n",rrr,resname,s[h],iii,x[h]/10,y[h]/10,z[h]/10
    }
  printf "%10.5f%10.5f%10.5f\n",a,b,c
      }' random_pos.dat > start_random.gro

