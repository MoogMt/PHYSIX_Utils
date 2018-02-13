program hydrogen
implicit none
integer, parameter :: np = 300 ! number of radial points
real(8), parameter :: a0 = 1._8 ! Borh's radius
real(8), dimension(:,:), allocatable :: rad
real(8), dimension(:), allocatable :: c
integer :: i, j, l, m, n, fact
real(8) :: r, dr, norm, theta, phi
real(8) :: laguerre, s_harm
character(9) :: flnm="hydro.0_0"

! Radial part
print '("Enter n: ",$)' ; read *, n
write(flnm(7:7),'(i1)') n
allocate(rad(0:n-1,0:np), c(0:n-1))
do l = 0, n-1
  c(l) =  dble(fact(n-l-1))/(2*n*fact(n+l)**3)*(2/(n*a0))**(l+1.5_8)
enddo
dr = 30._8/np
do i = 0, np
  r = i*dr
  do l = 0, n-1
    rad(l,i) =  c(l)*r**l*exp(-r/(n*a0))*laguerre(n+l-1,2*l+1,2*r/(n*a0))
  enddo
enddo
do l = 0, n-1
  norm = 0._8
  do i = 0, np
    norm = norm + (i*dr*rad(l,i))**2*dr
  enddo
  rad(l,:) = rad(l,:)/sqrt(norm)
enddo

! Write out with angular stuff
print '("Enter l: ",$)' ; read *, l
write(flnm(9:9),'(i1)') l
print '("Output file: ",a)', flnm
open(10,file=flnm)
phi = 0.
do j = 0, 100
 theta = j*acos(-1._8)/100
  do i = 0, np/2
    write(10,'(10g15.7)') i*dr, theta,  &
                          (rad(l,i)*s_harm(l,m,theta,phi), m = 0, l)
  enddo
  write(10,*)
enddo
close(10)
end program hydrogen

include 'laguerre.f90'
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

