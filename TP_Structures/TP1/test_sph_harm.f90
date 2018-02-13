program t_s_h
implicit none
real(8), parameter :: pi = acos(-1._8)
real(8) :: theta, phi, s_harm
integer :: l, m, ithe, iphi, fact, n
character(10) :: flnm= 'sph_harm_0'
print '("Enter l: ",$)' ; read *, l
write(flnm(10:10),'(i1)'), l
n = 45
open(10, file=flnm)
do ithe = 0, n
 theta = ithe*pi/n
 do iphi = 0, 2*n
 phi = iphi*pi/n
 write(10,'(10g13.5)') theta, phi, (s_harm(l,m,theta,phi), m = 0, l )
 enddo
 write(10,*)
enddo
close(10)
end program t_s_h


include 'sph_harm.f90'

integer function fact(l)
implicit none
integer :: k, l, i
i = 1
if ( l > 1 ) then
 do k = 2, l
  i = i*k
 enddo
endif
fact = i
end function fact

