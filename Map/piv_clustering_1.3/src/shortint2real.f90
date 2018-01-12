!------------------------------------------------------------------------------------------!
!                                                                                          !
!------------------------------------------------------------------------------------------!
function shortint2real(shortint,rmax,rmin)
  implicit none
  double precision::shortint2real,rmax,rmin,intsize
#ifdef FOURBYTES  
  double precision::shortintmin,shortintmax
  integer(kind=4)::real2shortint,shortint
  parameter(shortintmin=-2147483648.d0,shortintmax=2147483647.d0,intsize=4294967295.d0)
#else
  integer::shortintmin,shortintmax
  integer(kind=2)::real2shortint,shortint
  parameter(shortintmin=-32768,shortintmax=32767,intsize=65535.d0)
#endif
  shortint2real=(shortint-shortintmin)*(rmax-rmin)/intsize+rmin

  return
end function shortint2real
