real(8) function s_harm(l,m,theta,phi)
implicit none
integer :: l, m, fact
real(8) :: theta, phi, assoc_legendre
real(8), parameter :: pi = acos(-1._8)
if ( m==0 ) then
  s_harm = sqrt((2*l+1)/(4*pi))*assoc_legendre(l,0,cos(theta))
else
  s_harm = sqrt(2*(2*l+1)*fact(l-m)/(4*pi*fact(l+m)))* &
              assoc_legendre(l,m,cos(theta))*cos(m*phi)
endif
end function s_harm

include 'legendre.f90'
