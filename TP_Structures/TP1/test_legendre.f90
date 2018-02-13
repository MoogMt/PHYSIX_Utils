program test_leg
implicit none
integer :: k, l, m
real(8) :: x, dx, legendre, assoc_legendre
character(10) :: flnm = 'legendre_0'

dx = 0.01
print '("l=",$)' ; read *, l
write (flnm(10:10),'(i1)') l
print '("Output file: ",a)', flnm
open(10,file=flnm)
do k = -100, 100
x = k*dx
write(10,'(10g13.5)') x, legendre(l,x), (assoc_legendre(l,m,x), m = 0, l )
enddo
close(10)
end program test_leg

include 'legendre.f90'
