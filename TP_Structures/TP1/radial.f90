program radial
implicit none
integer, parameter :: np = 300
real(8), parameter :: a0 = 1._8
real(8), dimension(:,:), allocatable :: rad
real(8), dimension(:), allocatable :: c
integer :: i, n, l, fact
real(8) :: r, dr, norm
real(8) :: laguerre
character(10) :: flnm="radial.000"

print '("Enter n: ",$)' ; read *, n
write(flnm(8:10),'(i3.3)') n
print '("Output file: ",a)', flnm
open(10,file=flnm)
allocate(rad(0:n-1,0:np), c(0:n-1))
do l = 0, n-1
  c(l) =  dble(fact(n-l-1))/(2*n*fact(n+l)**3)*(2/(n*a0))**(l+1.5_8)
enddo
dr = 30._8/np
do i = 0, np
  r = i*dr
  do l = 0, n-1
    rad(l,i) =  c(l)*r**l*exp(-r/(n*a0))*laguerre(n-l-1,2*l+1,2*r/(n*a0))
  enddo
enddo
do l = 0, n-1
  norm = 0._8
  do i = 0, np
    norm = norm + (i*dr*rad(l,i))**2*dr
  enddo
  rad(l,:) = rad(l,:)/sqrt(norm)
enddo

do i = 0, np/2
  write(10,'(10g13.5)') i*dr, (rad(l,i), l = 0, n-1)
enddo
close(10)
end program radial

include 'laguerre.f90'


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

