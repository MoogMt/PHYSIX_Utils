!------------------------------------------------------------------------------------------!
!                                                                                          !
!------------------------------------------------------------------------------------------!
function real2shortint(r,rmax,rmin)
  implicit none
  double precision::r,rmax,rmin
#ifdef FOURBYTES  
  double precision::shortintmin,shortintmax,intsize
  integer(kind=4)::real2shortint
  parameter(shortintmin=-2147483648.d0,shortintmax=2147483647.d0,intsize=4294967295.d0)
#else
  integer::shortintmin,shortintmax,intsize
  integer(kind=2)::real2shortint
  parameter(shortintmin=-32768,shortintmax=32767,intsize=65535)
#endif
  real2shortint=nint((r-rmin)/(rmax-rmin)*(intsize)+shortintmin)
  return
end function real2shortint
!------------------------------------------------------------------------------------------!
